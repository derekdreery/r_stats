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
//! R is released as GPLv2 which I interpret as meaning this library must also be released as
//! GPLv2. If all the functions were replaced with native rust ones, then the license could be
//! changed to something more permissive.

/// This module provides the output of bindgen, in case the raw functions are useful.
pub mod ffi {
    include!(concat!(env!("OUT_DIR"), "/bindings.rs"));
}

// Normal distribution

/// Evaluate the probability density function of the normal distribution with mean `mu` and
/// variance `sigma`<sup>2</sup> at `x`.
///
/// The formula for this function is
/// <math>
///   <mfrac>
///     <mrow>
///       <mn>1</mn>
///     </mrow>
///     <mrow>
///       <mi>σ</mi>
///       <msqrt>
///         <mn>2π</mn>
///       </msqrt>
///     </mrow>
///   </mfrac>
///   <msup>
///     <mn>e</mn>
///     <mrow>
///       <mo>-</mo>
///       <mfrac>
///         <mn>1</mn>
///         <mn>2</mn>
///       </mfrac>
///       <mrow>
///         <mo>(</mo>
///         <mfrac>
///           <mrow>
///             <mi>x</mi>
///             <mo>-</mo>
///             <mn>μ</mn>
///           </mrow>
///           <mn>σ</mn>
///         </mfrac>
///         <mo>)</mo>
///       </row>
///     </mrow>
///   </msup>
/// </math>
/// where μ is the mean, and σ<sup>2</sup> is the variance.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn normal_pdf(x: f64, mu: f64, sigma: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnorm4(x, mu, sigma, c_bool(give_log)) }
}

/// Evaluate the culmulative density function of the normal distribution with mean `mu` and
/// variance `sigma` squared at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn normal_cdf(x: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pnorm5(x, mu, sigma, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the normal distribution with mean `mu` and
/// variance `sigma` squared at probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn normal_quantile(p: f64, mu: f64, sigma: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qnorm5(p, mu, sigma, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the culmulative density function of the normal distribution with mean `mu` and
/// variance `sigma` at `x`. Both integrals (`(-∞, x)` and `(x, ∞)`) are returned in that order.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn normal_cdf_both(x: f64, mu: f64, sigma: f64, log_p: bool) -> (f64, f64) {
    // Since we are using `mu` and `sigma` rather than 0 and 1, we have to check those values for
    // edge cases.
    if x.is_nan() || mu.is_nan() || sigma.is_nan() {
        let nan = x + mu + sigma;
        (nan, nan)
    } else if !x.is_finite() && mu == x {
        (f64::NAN, f64::NAN)
    } else if sigma <= 0.0 {
        if sigma < 0.0 {
            (f64::NAN, f64::NAN)
        } else {
            (zero(log_p), one(log_p))
        }
    } else {
        let x_canon = (x - mu) / sigma;
        if !x_canon.is_finite() {
            return if x < mu {
                (zero(log_p), one(log_p))
            } else {
                (one(log_p), zero(log_p))
            };
        }
        let x = x_canon;
        let mut lower: f64 = 0.0;
        let mut upper: f64 = 0.0;
        unsafe { ffi::pnorm_both(x, &mut lower, &mut upper, 2, c_bool(log_p)) }
        (lower, upper)
    }
}

// todo uniform distribution (it kinda feels too trivial to include).

// Gamma distribution

/// Evaluate the probability density function of the gamma distribution with given `shape` and
/// `scale` at `x`.
///
/// The shape is sometimes labelled `alpha`, and the scale is sometimes parameterised as
/// (`1 / lambda`).
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn gamma_pdf(x: f64, shape: f64, scale: f64, give_log: bool) -> f64 {
    unsafe { ffi::dgamma(x, shape, scale, c_bool(give_log)) }
}

/// Evaluate the culmulative density function of the gamma distribution with `shape` and `scale` at
/// `x`.
///
/// The shape is sometimes labelled `alpha`, and the scale is sometimes parameterised as
/// (`1 / lambda`).
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn gamma_cdf(x: f64, shape: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pgamma(x, shape, scale, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the gamma distribution with `shape` and `scale` at
/// probability `p`.
///
/// The shape is sometimes labelled `alpha`, and the scale is sometimes parameterised as
/// (`1 / lambda`).
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn gamma_quantile(p: f64, shape: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qgamma(p, shape, scale, c_bool(lower_tail), c_bool(log_p)) }
}

// Beta distribution

/// Evaluate the probability density function of the beta distribution with parameters `a` and `b`
/// at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn beta_pdf(x: f64, a: f64, b: f64, give_log: bool) -> f64 {
    unsafe { ffi::dbeta(x, a, b, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the beta distribution with parameters `a` and
/// `b` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn beta_cdf(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pbeta(x, a, b, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the beta distribution with parameters `a` and `b` at
/// probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn beta_quantile(p: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qbeta(p, a, b, c_bool(lower_tail), c_bool(log_p)) }
}

// Log-normal distribution

/// Evaluate the probability density function of the log-normal distribution with parameters
/// `mean_log` and `sd_log` at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn lognormal_pdf(x: f64, mean_log: f64, sd_log: f64, give_log: bool) -> f64 {
    unsafe { ffi::dlnorm(x, mean_log, sd_log, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the log-normal distribution with parameters
/// `mean_log` and `sd_log` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn lognormal_cdf(x: f64, mean_log: f64, sd_log: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::plnorm(x, mean_log, sd_log, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the log-normal distribution with parameters `mean_log` and
/// `sd_log` at probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn lognormal_quantile(
    p: f64,
    mean_log: f64,
    sd_log: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::qlnorm(p, mean_log, sd_log, c_bool(lower_tail), c_bool(log_p)) }
}

// Chi-squared distribution

/// Evaluate the probability density function of the chi-squared distribution with `df` degrees of
/// freedom at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn chi_squared_pdf(x: f64, df: f64, give_log: bool) -> f64 {
    unsafe { ffi::dchisq(x, df, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the chi-squared distribution with `df`
/// degrees of freedom at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn chi_squared_cdf(x: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pchisq(x, df, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the chi-squared distribution with `df` degrees of freedom at
/// probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn chi_squared_quantile(p: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qchisq(p, df, c_bool(lower_tail), c_bool(log_p)) }
}

// Non-central chi-squared distribution

/// Evaluate the probability density function of the non-central chi-squared distribution with `df`
/// degrees of freedom and non-centrality parameter `ncp` at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn noncentral_chi_squared_pdf(x: f64, df: f64, ncp: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnchisq(x, df, ncp, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the non-central chi-squared distribution with
/// `df` degrees of freedom and non-centrality parameter `ncp` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn noncentral_chi_squared_cdf(x: f64, df: f64, ncp: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pnchisq(x, df, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the non-central chi-squared distribution with `df` degrees of
/// freedom and non-centrality parameter `ncp` at probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn noncentral_chi_squared_quantile(
    p: f64,
    df: f64,
    ncp: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::qnchisq(p, df, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

// F distribution

/// Evaluate the probability density function of the f distribution with parameters `df1` and `df2`
/// at `x`. TODO I think df1 is the numerator degrees of freedom (when viewed as the ratio of two
/// chi-squared dists. Check and doc this.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn f_pdf(x: f64, df1: f64, df2: f64, give_log: bool) -> f64 {
    unsafe { ffi::df(x, df1, df2, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the f distribution with parameters `df1` and
/// `df2` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn f_cdf(x: f64, df1: f64, df2: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pf(x, df1, df2, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the f distribution with parameters `df1` and `df2` at
/// probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn f_quantile(p: f64, df1: f64, df2: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qf(p, df1, df2, c_bool(lower_tail), c_bool(log_p)) }
}

// Student's t distribution

/// Evaluate the probability density function of the student's t distribution with degrees of
/// freedom `df` at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn students_t_pdf(x: f64, df: f64, give_log: bool) -> f64 {
    unsafe { ffi::dt(x, df, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the student's t distribution with degrees of
/// freedom `df` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn students_t_cdf(x: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pt(x, df, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the student's t distribution with degrees of
/// freedom `df` at probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn students_t_quantile(p: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qt(p, df, c_bool(lower_tail), c_bool(log_p)) }
}

// Binomial distribution

/// Evaluate the probability density function of the binomial distribution with `n` trials and
/// probability of success `p` at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn binomial_pdf(x: f64, n: f64, p: f64, give_log: bool) -> f64 {
    unsafe { ffi::dbinom(x, n, p, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the binomial distribution with `n` trials and
/// probability of success `p` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn binomial_cdf(x: f64, n: f64, p: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pbinom(x, n, p, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the binomial distribution with `n` trials and
/// probability of success `pr` at probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn binomial_quantile(p: f64, n: f64, pr: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qbinom(p, n, pr, c_bool(lower_tail), c_bool(log_p)) }
}

// Ignoring multnomial because it only has random generation (and is spelt wrong - I bet it isn't
// used much).

// Cauchy distribution

/// Evaluate the probability density function of the Cauchy distribution with parameters `location`
/// and `scale` at `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn cauchy_pdf(x: f64, location: f64, scale: f64, give_log: bool) -> f64 {
    unsafe { ffi::dcauchy(x, location, scale, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the Cauchy distribution with parameters
/// `location` and `scale` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn cauchy_cdf(x: f64, location: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pcauchy(x, location, scale, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the Cauchy distribution with parameters
/// `location` and `scale` at probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn cauchy_quantile(p: f64, location: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qcauchy(p, location, scale, c_bool(lower_tail), c_bool(log_p)) }
}

// Exponential distribution

/// Evaluate the probability density function of the exponential distribution with given `scale` at
/// `x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn exponential_pdf(x: f64, scale: f64, give_log: bool) -> f64 {
    unsafe { ffi::dexp(x, scale, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the exponential distribution with given
/// `scale` at `x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn exponential_cdf(x: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pexp(x, scale, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the exponential distribution with given `scale` at
/// probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn exponential_quantile(p: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qexp(p, scale, c_bool(lower_tail), c_bool(log_p)) }
}

// Geometric distribution

/// Evaluate the probability density function of the geometric distribution with trial probability
/// `p` at `x`. `P(X=x) = p(1-p)^x`.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn geometric_pdf(x: f64, p: f64, give_log: bool) -> f64 {
    unsafe { ffi::dgeom(x, p, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the geometric distribution with trial
/// probability `p` at `x`. `P(X=x) = p(1-p)^x`.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn geometric_cdf(x: f64, p: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pgeom(x, p, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the geometric distribution with trial probability `prob` at
/// probability `p`. `P(X=x) = p(1-p)^x`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn geometric_quantile(p: f64, prob: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qgeom(p, prob, c_bool(lower_tail), c_bool(log_p)) }
}

// Hypergeometric distribution

/// Evaluate the probability density function of the hypergeometric distribution. If `samples`
/// samples were taken at random from a collection of `succ` successes and `fail` failures, then
/// this function evaluates the chance of `x` successes being in the sample.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn hypergeometric_pdf(x: f64, succ: f64, fail: f64, samples: f64, give_log: bool) -> f64 {
    unsafe { ffi::dhyper(x, succ, fail, samples, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the hypergeometric distribution. If `samples`
/// samples were taken at random from a collection of `succ` successes and `fail` failures, then
/// this function evaluates the chance of `≤x` successes being in the sample.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn hypergeometric_cdf(
    x: f64,
    succ: f64,
    fail: f64,
    samples: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::phyper(x, succ, fail, samples, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the hypergeometric distribution. If `samples`
/// samples were taken at random from a collection of `succ` successes and `fail` failures, then
/// this function evaluates the number of samples where there is `p` probability there are <= that
/// many successes.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn hypergeometric_quantile(
    p: f64,
    succ: f64,
    fail: f64,
    samples: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::qhyper(p, succ, fail, samples, c_bool(lower_tail), c_bool(log_p)) }
}

// Negative binomial distribution

/// Evaluate the probability density function of the negative binomial distribution. If there is
/// `prob` probability of successin a Bernoulli trial, and we perform trials until we see `size`
/// successes, then this function returns the probability of seeing `x` failures in those trials.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn neg_binomial_pdf(x: f64, size: f64, prob: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnbinom(x, size, prob, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the negative binomial distribution. If there
/// is `prob` probability of successin a Bernoulli trial, and we perform trials until we see `size`
/// successes, then this function returns the probability of seeing `≤x` failures in those trials.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn neg_binomial_cdf(x: f64, size: f64, prob: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pnbinom(x, size, prob, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the negative binomial distribution. If there
/// is `prob` probability of successin a Bernoulli trial, and we perform trials until we see `size`
/// successes, then this function returns the number of failures we would expect to see with
/// probability `p`.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn neg_binomial_quantile(p: f64, size: f64, prob: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qnbinom(p, size, prob, c_bool(lower_tail), c_bool(log_p)) }
}

// TODO negative binomial `_mu` - I don't know what it means or what it evaluates.

// Poisson distribution

/// Evaluate the probability density function of the poisson distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn poisson_pdf(x: f64, lambda: f64, give_log: bool) -> f64 {
    unsafe { ffi::dpois(x, lambda, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the poisson distribution.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn poisson_cdf(x: f64, lambda: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::ppois(x, lambda, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the poisson distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn poisson_quantile(p: f64, lambda: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qpois(p, lambda, c_bool(lower_tail), c_bool(log_p)) }
}

// Weibull distribution

/// Evaluate the probability density function of the Weibull distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn weibull_pdf(x: f64, shape: f64, scale: f64, give_log: bool) -> f64 {
    unsafe { ffi::dweibull(x, shape, scale, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the Weibull distribution.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn weibull_cdf(x: f64, shape: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pweibull(x, shape, scale, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the Weibull distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn weibull_quantile(p: f64, shape: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qweibull(p, shape, scale, c_bool(lower_tail), c_bool(log_p)) }
}

// Logistic distribution

/// Evaluate the probability density function of the logistic distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn logistic_pdf(x: f64, location: f64, scale: f64, give_log: bool) -> f64 {
    unsafe { ffi::dlogis(x, location, scale, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the logistic distribution.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn logistic_cdf(x: f64, location: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::plogis(x, location, scale, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the logistic distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn logistic_quantile(p: f64, location: f64, scale: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qlogis(p, location, scale, c_bool(lower_tail), c_bool(log_p)) }
}

// Non-central beta distribution

/// Evaluate the probability density function of the non-central beta distribution.
///
/// These functions use different algorithms to the beta distribution, so setting `ncp = 0` might
/// not give the same result.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn non_central_beta_pdf(x: f64, a: f64, b: f64, ncp: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnbeta(x, a, b, ncp, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the non-central beta distribution.
///
/// These functions use different algorithms to the beta distribution, so setting `ncp = 0` might
/// not give the same result.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn non_central_beta_cdf(
    x: f64,
    a: f64,
    b: f64,
    ncp: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::pnbeta(x, a, b, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the non-central beta distribution.
///
/// These functions use different algorithms to the beta distribution, so setting `ncp = 0` might
/// not give the same result.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn non_central_beta_quantile(
    p: f64,
    a: f64,
    b: f64,
    ncp: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::qnbeta(p, a, b, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

// Non-central f distribution

/// Evaluate the probability density function of the non-central f distribution.
///
/// These functions use different algorithms to the f distribution, so setting `ncp = 0` might
/// not give the same result.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn non_central_f_pdf(x: f64, df1: f64, df2: f64, ncp: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnf(x, df1, df2, ncp, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the non-central f distribution.
///
/// These functions use different algorithms to the f distribution, so setting `ncp = 0` might
/// not give the same result.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn non_central_f_cdf(
    x: f64,
    df1: f64,
    df2: f64,
    ncp: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::pnf(x, df1, df2, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the non-central f distribution.
///
/// These functions use different algorithms to the f distribution, so setting `ncp = 0` might
/// not give the same result.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn non_central_f_quantile(
    p: f64,
    df1: f64,
    df2: f64,
    ncp: f64,
    lower_tail: bool,
    log_p: bool,
) -> f64 {
    unsafe { ffi::qnf(p, df1, df2, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

// Non-central student's t distribution

/// Evaluate the probability density function of the non-central student's t distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn non_central_t_pdf(x: f64, df: f64, ncp: f64, give_log: bool) -> f64 {
    unsafe { ffi::dnt(x, df, ncp, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the non-central student's t distribution.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn non_central_t_cdf(x: f64, df: f64, ncp: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pnt(x, df, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the non-central student's t distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn non_central_t_quantile(p: f64, df: f64, ncp: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qnt(p, df, ncp, c_bool(lower_tail), c_bool(log_p)) }
}

// Studentized range distribution

/// Evaluate the culmulative distribution function of the studentized range distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn tukey_pdf(q: f64, rr: f64, cc: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::ptukey(q, rr, cc, df, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the studentized range distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn tukey_quantile(p: f64, rr: f64, cc: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qtukey(p, rr, cc, df, c_bool(lower_tail), c_bool(log_p)) }
}

// Wilcoxon rank-sum distribution

/// Evaluate the probability density function of the Wilcoxon rank-sum distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn wilcox_pdf(x: f64, m: f64, n: f64, give_log: bool) -> f64 {
    unsafe { ffi::dwilcox(x, m, n, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the Wilcoxon rank-sum distribution.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn wilcox_cdf(x: f64, m: f64, n: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::pwilcox(x, m, n, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the Wilcoxon rank-sum distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn wilcox_quantile(p: f64, m: f64, n: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qwilcox(p, m, n, c_bool(lower_tail), c_bool(log_p)) }
}

// Wilcoxon signed rank distribution

/// Evaluate the probability density function of the Wilcoxon signed rank distribution.
///
/// If `give_log` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn signrank_pdf(x: f64, n: f64, give_log: bool) -> f64 {
    unsafe { ffi::dsignrank(x, n, c_bool(give_log)) }
}

/// Evaluate the culmulative distribution function of the signrankon signed rank distribution.
///
/// If `lower_tail` is true, the integral from `-∞` to `x` is evaluated, else the
/// integral from `x` to `∞` is evaluated instead. "Usual" behaviour corresponds to
/// `true`. Using `lower_tail = false` gives higher numerical accuracy than performing the
/// calculation `1 - result` on `lower_tail = true` (when the result is close to 1).
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()` on the result.
pub fn signrank_cdf(x: f64, n: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::psignrank(x, n, c_bool(lower_tail), c_bool(log_p)) }
}

/// Evaluate the quantile function of the signrankon signed rank distribution.
///
/// If `lower_tail` is true, then `p` is the integral from `-∞` to `x`, else it is the integral
/// from `x` to `∞`. "Usual" behaviour corresponds to `true`.
///
/// If `log_p` is true, the natural logarithm of the value will be returned, with potentially
/// higher numerical accuracy than calling `.ln()`.
pub fn signrank_quantile(p: f64, n: f64, lower_tail: bool, log_p: bool) -> f64 {
    unsafe { ffi::qsignrank(p, n, c_bool(lower_tail), c_bool(log_p)) }
}

// TODO maybe r* functions (you can get their behaviour from the `rand` family
// of crates though). Also some general math functions like the beta/gamma functions, and others.
// I'm not sure if they should be included in this crate or not.

/// Helper to convert rust bools to c bools. Should be compiled away during optimization.
#[inline(always)]
fn c_bool(v: bool) -> i32 {
    if v {
        1
    } else {
        0
    }
}

/// Helper to return zero, taking logs if necessary
#[inline(always)]
fn zero(log: bool) -> f64 {
    if log {
        f64::NEG_INFINITY
    } else {
        0.0
    }
}

/// Helper to return 1, taking logs if necessary
#[inline(always)]
fn one(log: bool) -> f64 {
    if log {
        0.0
    } else {
        1.0
    }
}

#[cfg(test)]
mod tests {
    use super::ffi;
    use float_cmp::approx_eq;

    #[test]
    fn qf() {
        unsafe {
            // A not-very-accurate test, but it is from the R docs.
            assert!((ffi::qf(0.95, 5., 2., 1, 0) - 19.296).abs() < 1e-2);
        }
    }
}
