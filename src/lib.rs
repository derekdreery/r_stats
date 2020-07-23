#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
//! Things I've done to make R's nmath library work
//!   1. Copied the nmath source
//!   2. Set `MATHLIB_STANDALONE`: telling nmath to build for use outside R
//!   3. Copied R's includes since they are still needed even with `MATHLIB_STANDALONE`.
//!   4. Copy `Rconfig.h` and `config.h` from the mingw dir (normally these files would be
//!      generated, but I think the math libs use very little of them).
//!   5. Copy `Rmath.h` from `Rmath.h0.in` and replace the build system placeholders (there are
//!      only 2).
//!   6. Manually whitelist the functions to include, otherwise we pull in a lot of stuff,
//!      including stuff outside the source tree (and also it doesn't compile for some reason).
//!   7. `Rmath.h` wierdly sets the normal functions as aliases - we have to use the raw function
//!      name here e.g. `dnorm4`.
//!
//! I've done some work on rust versions of these functions in the `riir` branch. This branch
//! should be used for testing those functions in future.
//!
//! R is released as GPLv2 which I interperet as meaning this library must also be released as
//! GPLv2. If all the functions were replaced with native rust ones, then the license could be
//! changed to something more permissive.

/// This module provides the output of bindgen, in case the raw functions are useful.
pub mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

/// Evaluate the probability density function of the normal distribution with mean `mu` and
/// variance `sigma` at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than naively calling `.ln()`.
pub fn normal_pdf(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnorm4(x, mu, sigma, c_bool(give_log)) }
}

/// Evaluate the culmulative density function of the normal distribution with mean `mu` and
/// variance `sigma` at `x`. (The integral of the pdf up to `x`).
///
/// If `lower_tail` is true, the integral from `-f64::NEG_INFINITY` to `x` is evaluated, else the
/// integral from `x` to `f64::INFINITY` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Numerical accuracy may be higher than naively calculating `1.0 - normal_cdf(.., false,
/// ..)` for `normal_cdf(.., true, ..)`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than naively calling `.ln()`.
pub fn normal_cdf(x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pnorm5(x, mu, sigma, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the normal distribution with mean `mu` and
/// variance `sigma` at probability `p`. The quantile function is the inverse of the cdf.
///
/// If `lower_tail` is true, the integral from `-f64::NEG_INFINITY` to `x` is evaluated, else the
/// integral from `x` to `f64::INFINITY` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Numerical accuracy may be higher than naively calculating `1.0 - normal_cdf(.., false,
/// ..)` for `normal_cdf(.., true, ..)`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than naively calling `.ln()`.
pub fn normal_quantile(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qnorm5(p, mu, sigma, c_bool(lower_tail), c_bool(log_p)) }
}

// TODO all other distributions.

/// Helper to convert rust bools to c bools. Should be compiled away during optimization.
#[inline(always)]
fn c_bool(v: bool) -> i32 {
    if v {
        1
    } else {
        0
    }
}

#[cfg(test)]
mod tests {
    use super::ffi;
    use float_cmp::approx_eq;

    #[test]
    fn df() {
        unsafe {
            // A not-very-accurate test, but it is from the R docs.
            assert!((ffi::qf(0.95, 5., 2., 1, 0) - 19.296).abs() < 1e-2);
        }
    }
}
