! rustprop from Fortran, via iso_c_binding.
!
! The C ABI is the interoperable one, so no shim library is needed — Fortran
! calls it directly.
!
!   gfortran rustprop.f90 -o rustprop_demo -L<sdk>/lib -lrustprop
!   LD_LIBRARY_PATH=<sdk>/lib ./rustprop_demo

module rustprop
    use, intrinsic :: iso_c_binding
    implicit none

    integer(c_int), parameter :: RUSTPROP_OK = 0
    integer(c_int), parameter :: RUSTPROP_UNAVAILABLE = 102

    interface
        integer(c_int) function c_props_si(output, n1, v1, n2, v2, fluid, out) &
                bind(C, name="rustprop_props_si")
            import :: c_char, c_double, c_int
            character(kind=c_char), dimension(*), intent(in) :: output, n1, n2, fluid
            real(c_double), value :: v1, v2
            real(c_double), intent(out) :: out
        end function

        integer(c_int) function c_ha_props_si(output, n1, v1, n2, v2, n3, v3, out) &
                bind(C, name="rustprop_ha_props_si")
            import :: c_char, c_double, c_int
            character(kind=c_char), dimension(*), intent(in) :: output, n1, n2, n3
            real(c_double), value :: v1, v2, v3
            real(c_double), intent(out) :: out
        end function

        integer(c_int) function c_props_si_many(output, n1, v1, n2, v2, fluid, n, out) &
                bind(C, name="rustprop_props_si_many")
            import :: c_char, c_double, c_int, c_size_t
            character(kind=c_char), dimension(*), intent(in) :: output, n1, n2, fluid
            real(c_double), dimension(*), intent(in) :: v1, v2
            integer(c_size_t), value :: n
            real(c_double), dimension(*), intent(out) :: out
        end function

        integer(c_size_t) function c_last_error_message(buf, len) &
                bind(C, name="rustprop_last_error_message")
            import :: c_char, c_size_t
            character(kind=c_char), dimension(*), intent(out) :: buf
            integer(c_size_t), value :: len
        end function

        type(c_ptr) function c_backends() bind(C, name="rustprop_backends")
            import :: c_ptr
        end function

        type(c_ptr) function c_version() bind(C, name="rustprop_version")
            import :: c_ptr
        end function

        type(c_ptr) function c_upstream_version() bind(C, name="rustprop_upstream_version")
            import :: c_ptr
        end function

        integer(c_int) function c_has_backend(name) bind(C, name="rustprop_has_backend")
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: name
        end function

        integer(c_size_t) function c_fluid_count() bind(C, name="rustprop_fluid_count")
            import :: c_size_t
        end function
    end interface

contains

    !> A NUL-terminated copy of `s`, which is what the C side expects.
    function cstr(s) result(r)
        character(len=*), intent(in) :: s
        character(kind=c_char, len=:), allocatable :: r
        r = trim(s)//c_null_char
    end function

    !> Read a `const char *` the library returned. Those pointers are static,
    !> so nothing here owns or frees anything.
    function from_cptr(p) result(r)
        type(c_ptr), intent(in) :: p
        character(len=:), allocatable :: r
        character(kind=c_char), pointer :: chars(:)
        integer :: i, n
        if (.not. c_associated(p)) then
            r = ""
            return
        end if
        call c_f_pointer(p, chars, [huge(0_c_int)])
        n = 0
        do i = 1, huge(0_c_int)
            if (chars(i) == c_null_char) exit
            n = n + 1
        end do
        allocate (character(len=n) :: r)
        do i = 1, n
            r(i:i) = chars(i)
        end do
    end function

    !> The message for whatever just failed on this thread.
    function last_error() result(r)
        character(len=:), allocatable :: r
        character(kind=c_char), allocatable :: buf(:)
        integer(c_size_t) :: need, ignored
        integer :: i
        need = c_last_error_message(buf, 0_c_size_t)
        allocate (buf(need + 1))
        ignored = c_last_error_message(buf, need + 1)
        allocate (character(len=int(need)) :: r)
        do i = 1, int(need)
            r(i:i) = buf(i)
        end do
    end function

    !> PropsSI. `status` is RUSTPROP_OK on success; `value` is untouched otherwise.
    subroutine props_si(output, n1, v1, n2, v2, fluid, value, status)
        character(len=*), intent(in) :: output, n1, n2, fluid
        real(c_double), intent(in) :: v1, v2
        real(c_double), intent(out) :: value
        integer(c_int), intent(out) :: status
        status = c_props_si(cstr(output), cstr(n1), v1, cstr(n2), v2, cstr(fluid), value)
    end subroutine

    subroutine ha_props_si(output, n1, v1, n2, v2, n3, v3, value, status)
        character(len=*), intent(in) :: output, n1, n2, n3
        real(c_double), intent(in) :: v1, v2, v3
        real(c_double), intent(out) :: value
        integer(c_int), intent(out) :: status
        status = c_ha_props_si(cstr(output), cstr(n1), v1, cstr(n2), v2, &
                               cstr(n3), v3, value)
    end subroutine

    !> One output over many states; a failing state comes back as NaN.
    subroutine props_si_many(output, n1, v1, n2, v2, fluid, values, status)
        character(len=*), intent(in) :: output, n1, n2, fluid
        real(c_double), intent(in) :: v1(:), v2(:)
        real(c_double), intent(out) :: values(:)
        integer(c_int), intent(out) :: status
        status = c_props_si_many(cstr(output), cstr(n1), v1, cstr(n2), v2, &
                                 cstr(fluid), int(size(v1), c_size_t), values)
    end subroutine

    logical function has_backend(name)
        character(len=*), intent(in) :: name
        has_backend = (c_has_backend(cstr(name)) == 1)
    end function

end module rustprop

program rustprop_demo
    use rustprop
    use, intrinsic :: iso_c_binding
    implicit none

    real(c_double) :: d, h, w, one, t(3), p(3), many(3)
    integer(c_int) :: status
    integer :: failures, i

    failures = 0

    print '(a)', "rustprop "//from_cptr(c_version())// &
        " (CoolProp "//from_cptr(c_upstream_version())//")"
    print '(a)', "backends: "//from_cptr(c_backends())
    print '(a,i0)', "fluids compiled in: ", c_fluid_count()

    if (has_backend("heos")) then
        call props_si("Dmolar", "T", 300.0_c_double, "P", 101325.0_c_double, "Water", &
                      d, status)
        call check(status == RUSTPROP_OK .and. abs((d - 55317.35277350119_c_double)/d) &
                   < 1.0e-8_c_double, "PropsSI Dmolar Water")

        t = [300.0_c_double, 400.0_c_double, 500.0_c_double]
        p = [101325.0_c_double, 101325.0_c_double, 101325.0_c_double]
        call props_si_many("Dmolar", "T", t, "P", p, "Water", many, status)
        block
            logical :: same
            same = (status == RUSTPROP_OK)
            do i = 1, 3
                call props_si("Dmolar", "T", t(i), "P", p(i), "Water", one, status)
                if (many(i) /= one) same = .false.
            end do
            call check(same, "batch equals scalar exactly")
        end block
    end if

    if (has_backend("if97")) then
        call props_si("H", "T", 300.0_c_double, "P", 101325.0_c_double, "IF97::Water", &
                      h, status)
        call check(status == RUSTPROP_OK .and. abs((h - 112665.04341853978_c_double)/h) &
                   < 1.0e-11_c_double, "PropsSI H IF97::Water")
    end if

    if (has_backend("humid-air")) then
        call ha_props_si("W", "T", 300.0_c_double, "P", 101325.0_c_double, &
                         "R", 0.5_c_double, w, status)
        call check(status == RUSTPROP_OK .and. w > 0.0_c_double .and. &
                   w < 1.0_c_double, "HAPropsSI W")
    end if

    call props_si("Dmolar", "T", 300.0_c_double, "P", 101325.0_c_double, "NoSuchFluid", &
                  d, status)
    call check(status /= RUSTPROP_OK .and. index(last_error(), "NoSuchFluid") > 0, &
               "a bad fluid fails, naming the key")

    print '(a)', ""
    if (failures == 0) then
        print '(a)', "PASSED (0 failures)"
    else
        print '(a,i0,a)', "FAILED (", failures, " failures)"
        stop 1
    end if

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        character(len=52) :: padded
        padded = what
        if (ok) then
            print '(a,a,a)', "  ", padded, " ok"
        else
            print '(a,a,a)', "  ", padded, " FAILED"
            failures = failures + 1
        end if
    end subroutine

end program rustprop_demo
