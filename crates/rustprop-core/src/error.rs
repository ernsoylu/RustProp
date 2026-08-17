//! Error type mirroring upstream `include/CoolProp/Exceptions.h` @ v8.0.0.
//!
//! One variant per upstream `ErrCode` that a thermodynamic code path can
//! raise. The three host-integration codes (`eHandle`, `eUnableToLoad`,
//! `eDirectorySize` — C handle API, shared-library loading, REFPROP paths)
//! have no equivalent in this port and are deliberately omitted; add them
//! only if a ported code path actually needs them.
//!
//! Fidelity note (PLAN.md): error *conditions* match upstream; message
//! strings are diagnostics, not a compatibility surface. Like upstream's
//! `what()`, `Display` shows only the message.

/// Error condition, mirroring upstream `CoolPropBaseError::ErrCode`.
///
/// `#[non_exhaustive]`: semver headroom for tracking upstream (the omitted
/// host-integration codes, or codes a future CoolProp release adds). By
/// explicit ruling, `Param`/`InputPair`/`Phase` stay exhaustive — there,
/// compiler-forced match exhaustiveness is a port-completeness tool.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    /// `NotImplementedError`
    NotImplemented(String),
    /// `SolutionError` — an iterative solver failed to converge
    Solution(String),
    /// `AttributeError`
    Attribute(String),
    /// `OutOfRangeError` — input outside the range of validity
    OutOfRange(String),
    /// `ValueError` — invalid parameter/fluid/input names and values
    Value(String),
    /// `WrongFluidError`
    WrongFluid(String),
    /// `CompositionError`
    Composition(String),
    /// `InputError`
    Input(String),
    /// `NotAvailableError` — property not available for this fluid/model
    NotAvailable(String),
    /// `KeyError`
    Key(String),
    /// `MultipleSolutionsError` — a saturation flash input maps to more than
    /// one state (upstream GH #2773)
    MultipleSolutions(String),
}

impl Error {
    /// The diagnostic message (upstream `what()`).
    pub fn message(&self) -> &str {
        match self {
            Error::NotImplemented(m)
            | Error::Solution(m)
            | Error::Attribute(m)
            | Error::OutOfRange(m)
            | Error::Value(m)
            | Error::WrongFluid(m)
            | Error::Composition(m)
            | Error::Input(m)
            | Error::NotAvailable(m)
            | Error::Key(m)
            | Error::MultipleSolutions(m) => m,
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.message())
    }
}

impl std::error::Error for Error {}

pub type Result<T> = std::result::Result<T, Error>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn display_is_message_only_like_upstream_what() {
        let variants = [
            Error::NotImplemented("m".into()),
            Error::Solution("m".into()),
            Error::Attribute("m".into()),
            Error::OutOfRange("m".into()),
            Error::Value("m".into()),
            Error::WrongFluid("m".into()),
            Error::Composition("m".into()),
            Error::Input("m".into()),
            Error::NotAvailable("m".into()),
            Error::Key("m".into()),
            Error::MultipleSolutions("m".into()),
        ];
        for e in variants {
            assert_eq!(e.to_string(), "m");
            assert_eq!(e.message(), "m");
        }
    }
}
