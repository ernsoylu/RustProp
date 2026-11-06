//! A compact binary form of the fluid-data types, for loading fluids at run
//! time instead of linking them (PLAN.md Decisions, "Runtime fluids").
//!
//! The engines take `&'static` data, so a decoded fluid is leaked once and
//! lives for the process, exactly like a linked one. A registry must decode
//! each fluid once: the facade's caches are keyed by the data's address.
//!
//! Format (little-endian): an 8-byte magic naming the record kind and format
//! version, then the value. `f64` is its IEEE bit pattern, so NaN and the
//! infinities in the data survive exactly; `bool` and tags are one byte;
//! lengths and `usize` are `u32`/`u64`; strings are UTF-8. Trailing bytes are
//! refused.
//!
//! Exactness by construction: one macro writes both directions of each type
//! from one field list, in one order, and encoding destructures every value
//! without `..`, so a field cannot be left out of the format without a
//! compile error. With that, `encode(decode(encode(x))) == encode(x)` proves
//! the decoded value equal to `x` field for field, bit for bit; the golden
//! suite checks it for every fluid and the mixture tables.

use crate::fluid::*;
use crate::{Error, Result};

const FLUID_MAGIC: &[u8; 8] = b"RPFLUID1";
const MIXTURES_MAGIC: &[u8; 8] = b"RPMIXTB1";

/// The mixture tables `MixtureModel` and the predefined-mixture lookup read.
pub struct MixtureTables {
    pub binary_pairs: &'static [MixBinaryPair],
    pub departure_fns: &'static [MixDepartureFn],
    pub predefined: &'static [PredefinedMixture],
}

/// Encodes one fluid.
pub fn encode_fluid(data: &FluidData) -> Vec<u8> {
    let mut w = FLUID_MAGIC.to_vec();
    data.put(&mut w);
    w
}

/// Decodes one fluid and leaks it, so it is `&'static` like a linked one.
pub fn decode_fluid(bytes: &[u8]) -> Result<&'static FluidData> {
    let mut r = Reader::new(bytes, FLUID_MAGIC)?;
    let data = FluidData::take(&mut r)?;
    r.finish()?;
    Ok(Box::leak(Box::new(data)))
}

/// The name of the encoded fluid, read without decoding (or leaking) it.
pub fn fluid_name(bytes: &[u8]) -> Result<&str> {
    let mut r = Reader::new(bytes, FLUID_MAGIC)?;
    r.str()
}

/// Encodes the three mixture tables.
pub fn encode_mixtures(t: &MixtureTables) -> Vec<u8> {
    let mut w = MIXTURES_MAGIC.to_vec();
    t.binary_pairs.put(&mut w);
    t.departure_fns.put(&mut w);
    t.predefined.put(&mut w);
    w
}

/// Decodes the mixture tables and leaks them.
pub fn decode_mixtures(bytes: &[u8]) -> Result<&'static MixtureTables> {
    let mut r = Reader::new(bytes, MIXTURES_MAGIC)?;
    let t = MixtureTables {
        binary_pairs: Pack::take(&mut r)?,
        departure_fns: Pack::take(&mut r)?,
        predefined: Pack::take(&mut r)?,
    };
    r.finish()?;
    Ok(Box::leak(Box::new(t)))
}

fn bad(what: &str) -> Error {
    Error::Value(format!("packed fluid data is malformed: {what}"))
}

struct Reader<'a> {
    bytes: &'a [u8],
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8], magic: &[u8; 8]) -> Result<Self> {
        match bytes.split_at_checked(8) {
            Some((head, rest)) if head == magic => Ok(Reader { bytes: rest }),
            _ => Err(Error::Value(format!(
                "not packed fluid data of this kind and version (expected {})",
                String::from_utf8_lossy(magic)
            ))),
        }
    }
    fn take_bytes(&mut self, n: usize) -> Result<&'a [u8]> {
        let (head, rest) = self
            .bytes
            .split_at_checked(n)
            .ok_or_else(|| bad("truncated"))?;
        self.bytes = rest;
        Ok(head)
    }
    fn array<const N: usize>(&mut self) -> Result<[u8; N]> {
        Ok(self.take_bytes(N)?.try_into().expect("length checked"))
    }
    fn u8(&mut self) -> Result<u8> {
        Ok(self.array::<1>()?[0])
    }
    /// A length, refused when the rest of the input could not hold that many
    /// elements of at least `min` bytes each (so a corrupt length cannot ask
    /// for a huge allocation).
    fn len(&mut self, min: usize) -> Result<usize> {
        let n = u32::from_le_bytes(self.array()?) as usize;
        if n.saturating_mul(min) > self.bytes.len() {
            return Err(bad("length exceeds the data"));
        }
        Ok(n)
    }
    fn str(&mut self) -> Result<&'a str> {
        let n = self.len(1)?;
        std::str::from_utf8(self.take_bytes(n)?).map_err(|_| bad("string is not UTF-8"))
    }
    fn finish(self) -> Result<()> {
        if self.bytes.is_empty() {
            Ok(())
        } else {
            Err(bad("trailing bytes"))
        }
    }
}

fn put_len(w: &mut Vec<u8>, n: usize) {
    let n = u32::try_from(n).expect("a packed length fits in u32");
    w.extend_from_slice(&n.to_le_bytes());
}

trait Pack: Sized {
    /// The fewest bytes one encoded value takes.
    const MIN: usize;
    fn put(&self, w: &mut Vec<u8>);
    fn take(r: &mut Reader<'_>) -> Result<Self>;
}

impl Pack for f64 {
    const MIN: usize = 8;
    fn put(&self, w: &mut Vec<u8>) {
        w.extend_from_slice(&self.to_bits().to_le_bytes());
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        Ok(f64::from_bits(u64::from_le_bytes(r.array()?)))
    }
}

impl Pack for bool {
    const MIN: usize = 1;
    fn put(&self, w: &mut Vec<u8>) {
        w.push(u8::from(*self));
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        match r.u8()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(bad("bool is not 0 or 1")),
        }
    }
}

impl Pack for usize {
    const MIN: usize = 8;
    fn put(&self, w: &mut Vec<u8>) {
        w.extend_from_slice(&(*self as u64).to_le_bytes());
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        usize::try_from(u64::from_le_bytes(r.array()?)).map_err(|_| bad("usize overflow"))
    }
}

impl Pack for &'static str {
    const MIN: usize = 4;
    fn put(&self, w: &mut Vec<u8>) {
        put_len(w, self.len());
        w.extend_from_slice(self.as_bytes());
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        Ok(r.str()?.to_owned().leak())
    }
}

impl<T: Pack + 'static> Pack for &'static [T] {
    const MIN: usize = 4;
    fn put(&self, w: &mut Vec<u8>) {
        put_len(w, self.len());
        for x in self.iter() {
            x.put(w);
        }
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        let n = r.len(T::MIN)?;
        let mut v = Vec::with_capacity(n.min(r.bytes.len()));
        for _ in 0..n {
            v.push(T::take(r)?);
        }
        Ok(v.leak())
    }
}

impl<T: Pack> Pack for Option<T> {
    const MIN: usize = 1;
    fn put(&self, w: &mut Vec<u8>) {
        match self {
            None => w.push(0),
            Some(x) => {
                w.push(1);
                x.put(w);
            }
        }
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        match r.u8()? {
            0 => Ok(None),
            1 => Ok(Some(T::take(r)?)),
            _ => Err(bad("option tag")),
        }
    }
}

/// A struct: every field, in the listed order. The `let` without `..` is what
/// makes a forgotten field a compile error.
macro_rules! pack_struct {
    ($($ty:ident { $($f:ident),* $(,)? })*) => {$(
        impl Pack for $ty {
            // Every packed struct holds at least one non-empty field.
            const MIN: usize = 1;
            fn put(&self, w: &mut Vec<u8>) {
                let $ty { $($f),* } = self;
                $( $f.put(w); )*
            }
            fn take(r: &mut Reader<'_>) -> Result<Self> {
                Ok($ty { $($f: Pack::take(r)?),* })
            }
        }
    )*};
}

/// An enum of unit and struct variants: a one-byte tag, then the variant's
/// fields. The match has no wildcard, so a new variant is a compile error.
macro_rules! pack_enum {
    ($($ty:ident { $($tag:literal => $v:ident { $($f:ident),* $(,)? }),* $(,)? })*) => {$(
        impl Pack for $ty {
            const MIN: usize = 1;
            fn put(&self, w: &mut Vec<u8>) {
                match self {
                    $( $ty::$v { $($f),* } => {
                        w.push($tag);
                        $( $f.put(w); )*
                    } )*
                }
            }
            fn take(r: &mut Reader<'_>) -> Result<Self> {
                match r.u8()? {
                    $( $tag => Ok($ty::$v { $($f: Pack::take(r)?),* }), )*
                    _ => Err(bad(concat!(stringify!($ty), " tag"))),
                }
            }
        }
    )*};
}

pack_struct! {
    FluidData { name, cas, aliases, environmental, eos, ancillaries, states, transport }
    Environmental { gwp20, gwp100, gwp500, odp, hh, ph, fh }
    Transport { viscosity, conductivity }
    Conductivity { dilute, residual, critical }
    Viscosity { epsilon_over_k, sigma_eta, dilute, initial_density, higher_order }
    Eos {
        gas_constant, molar_mass, p_max, t_max, t_triple, acentric, pseudo_pure, reducing,
        sat_min_liquid, sat_min_vapor, hs_anchor, max_sat_t, max_sat_p, alpha0, alphar,
        superancillary,
    }
    SuperAncillaryData { p, rho_l, rho_v, t_crit_num, rho_crit_num, check_points }
    ChebyshevInterval { xmin, xmax, coef }
    SuperAncCheckPoint { t, p, rho_l, rho_v, p_ratio, rho_l_ratio, rho_v_ratio }
    StatePoint { t, p, rhomolar, hmolar, smolar }
    Ancillaries {
        p_s, p_v_split, rho_l, rho_v, h_l, h_lv, s_l, s_lv, surface_tension, melting_line,
    }
    RationalPolyAncillary { a, b, max_abs_error, t_min, t_max }
    MeltingLine { t_m, kind }
    SimonMeltPart { t_0, a, c, p_0, t_min, t_max }
    PolyMeltPart { t_0, p_0, t_min, t_max, a, t }
    SurfaceTension { a, n, tc }
    SaturationAncillary { anc_type, n, t, t_r, reducing_value, using_tau_r, t_min, t_max }
    States { critical, triple_liquid, triple_vapor }
    MixBinaryPair { cas1, cas2, beta_t, gamma_t, beta_v, gamma_v, f, function }
    PredefinedMixture { name, fluids, mole_fractions }
    MixDepartureFn { name, kind, npower, n, d, t, l, eta, epsilon, beta, gamma }
}

pack_enum! {
    ConductivityDilute {
        0 => RatioOfPolynomials { a, n, b, m, t_reducing },
        1 => Eta0AndPoly { a, t },
        2 => Hardcoded { name },
    }
    ConductivityResidual {
        0 => Polynomial { b, t, d, t_reducing, rhomass_reducing },
        1 => PolynomialAndExponential { a, t, d, gamma, l },
    }
    ConductivityCritical {
        0 => SimplifiedOlchowySengers { k, r0, gamma, nu, big_gamma, zeta0, qd, t_ref },
        1 => Hardcoded { name },
    }
    ViscosityDilute {
        0 => KineticTheory {},
        1 => CollisionIntegral { a, t, c, molar_mass },
        2 => PowersOfT { a, t },
        3 => PowersOfTr { a, t, t_reducing },
        4 => CollisionIntegralPowersOfTstar { a, t, c, t_reducing },
        5 => Hardcoded { name },
    }
    ViscosityInitialDensity {
        0 => RainwaterFriend { b, t },
        1 => Empirical { n, d, t, t_reducing, rhomolar_reducing },
    }
    ViscosityHigherOrder {
        0 => ModifiedBatschinskiHildebrand {
            a, d1, t1, gamma, l, f, d2, t2, g, h, p, q, t_reduce, rhomolar_reduce,
        },
        1 => FrictionTheory {
            ai, aa, ar, aaa, arr, adrdr, aii, arrr, aaaa, na, naa, nr, nrr, nii, nrrr, naaa,
            c1, c2, t_reduce,
        },
        2 => Hardcoded { name },
    }
    Alpha0Term {
        0 => Lead { a1, a2 },
        1 => LogTau { a },
        2 => PlanckEinstein { n, t },
        3 => PlanckEinsteinFunctionT { n, v, tcrit },
        4 => EnthalpyEntropyOffset { a1, a2, reference },
        5 => Power { n, t },
        6 => PlanckEinsteinGeneralized { n, t, c, d },
        7 => Cp0Constant { cp_over_r, tc, t0 },
        8 => Cp0PolyT { c, t, tc, t0 },
        9 => Cp0AlyLee { c, tc, t0 },
    }
    AlpharTerm {
        0 => Power { n, d, t, l },
        1 => Gaussian { n, d, t, eta, beta, gamma, epsilon },
        2 => NonAnalytic { n, a, b, beta, big_a, big_b, big_c, big_d },
        3 => Exponential { n, d, t, g, l },
        4 => DoubleExponential { n, d, t, gd, ld, gt, lt },
        5 => Lemmon2005 { n, d, t, l, m },
        6 => GaoB { n, t, d, eta, beta, gamma, epsilon, b },
    }
    MeltingLineKind {
        0 => Simon { parts },
        1 => PolynomialInTr { parts },
        2 => PolynomialInTheta { parts },
    }
    MixDepartureKind {
        0 => Gerg2008 {},
        1 => Exponential {},
        2 => GaussianExponential {},
    }
}

// The three enums with a tuple variant, which the macro does not take.

impl Pack for ViscosityModel {
    const MIN: usize = 1;
    fn put(&self, w: &mut Vec<u8>) {
        match self {
            ViscosityModel::Structured(v) => {
                w.push(0);
                v.put(w);
            }
            ViscosityModel::Hardcoded { name } => {
                w.push(1);
                name.put(w);
            }
            ViscosityModel::Chung {
                rhomolar_critical,
                acentric,
                molar_mass,
                t_critical,
                dipole_moment_d,
                kappa,
            } => {
                w.push(2);
                for x in [
                    rhomolar_critical,
                    acentric,
                    molar_mass,
                    t_critical,
                    dipole_moment_d,
                    kappa,
                ] {
                    x.put(w);
                }
            }
            ViscosityModel::Ecs {
                reference_fluid,
                psi_a,
                psi_t,
                psi_rhomolar_reducing,
                sigma_eta,
                epsilon_over_k,
            } => {
                w.push(3);
                reference_fluid.put(w);
                psi_a.put(w);
                psi_t.put(w);
                psi_rhomolar_reducing.put(w);
                sigma_eta.put(w);
                epsilon_over_k.put(w);
            }
            ViscosityModel::RhosrCs {
                c,
                c_liq,
                c_vap,
                rhosr_critical,
                x_crossover,
            } => {
                w.push(4);
                c.put(w);
                c_liq.put(w);
                c_vap.put(w);
                rhosr_critical.put(w);
                x_crossover.put(w);
            }
        }
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        Ok(match r.u8()? {
            0 => ViscosityModel::Structured(Pack::take(r)?),
            1 => ViscosityModel::Hardcoded {
                name: Pack::take(r)?,
            },
            2 => ViscosityModel::Chung {
                rhomolar_critical: Pack::take(r)?,
                acentric: Pack::take(r)?,
                molar_mass: Pack::take(r)?,
                t_critical: Pack::take(r)?,
                dipole_moment_d: Pack::take(r)?,
                kappa: Pack::take(r)?,
            },
            3 => ViscosityModel::Ecs {
                reference_fluid: Pack::take(r)?,
                psi_a: Pack::take(r)?,
                psi_t: Pack::take(r)?,
                psi_rhomolar_reducing: Pack::take(r)?,
                sigma_eta: Pack::take(r)?,
                epsilon_over_k: Pack::take(r)?,
            },
            4 => ViscosityModel::RhosrCs {
                c: Pack::take(r)?,
                c_liq: Pack::take(r)?,
                c_vap: Pack::take(r)?,
                rhosr_critical: Pack::take(r)?,
                x_crossover: Pack::take(r)?,
            },
            _ => return Err(bad("ViscosityModel tag")),
        })
    }
}

impl Pack for ConductivityModel {
    const MIN: usize = 1;
    fn put(&self, w: &mut Vec<u8>) {
        match self {
            ConductivityModel::Structured(c) => {
                w.push(0);
                c.put(w);
            }
            ConductivityModel::Hardcoded { name } => {
                w.push(1);
                name.put(w);
            }
            ConductivityModel::Ecs {
                reference_fluid,
                psi_a,
                psi_t,
                psi_rhomolar_reducing,
                f_int_a,
                f_int_t,
                f_int_t_reducing,
            } => {
                w.push(2);
                reference_fluid.put(w);
                psi_a.put(w);
                psi_t.put(w);
                psi_rhomolar_reducing.put(w);
                f_int_a.put(w);
                f_int_t.put(w);
                f_int_t_reducing.put(w);
            }
        }
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        Ok(match r.u8()? {
            0 => ConductivityModel::Structured(Pack::take(r)?),
            1 => ConductivityModel::Hardcoded {
                name: Pack::take(r)?,
            },
            2 => ConductivityModel::Ecs {
                reference_fluid: Pack::take(r)?,
                psi_a: Pack::take(r)?,
                psi_t: Pack::take(r)?,
                psi_rhomolar_reducing: Pack::take(r)?,
                f_int_a: Pack::take(r)?,
                f_int_t: Pack::take(r)?,
                f_int_t_reducing: Pack::take(r)?,
            },
            _ => return Err(bad("ConductivityModel tag")),
        })
    }
}

impl<T: Pack + 'static> Pack for TransportModel<T> {
    const MIN: usize = 1;
    fn put(&self, w: &mut Vec<u8>) {
        match self {
            TransportModel::Absent => w.push(0),
            TransportModel::Unported => w.push(1),
            TransportModel::Model(m) => {
                w.push(2);
                m.put(w);
            }
        }
    }
    fn take(r: &mut Reader<'_>) -> Result<Self> {
        Ok(match r.u8()? {
            0 => TransportModel::Absent,
            1 => TransportModel::Unported,
            2 => TransportModel::Model(T::take(r)?),
            _ => return Err(bad("TransportModel tag")),
        })
    }
}
