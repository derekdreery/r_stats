//! This module provides functions from R's `stats` module, but rewritten for rust.
//!
//! I'm sticking to the `stats` module's naming for functions, even though they aren't at all
//! clear. Maybe this could change in the future.
//!
//! I'm also not implementing the non-stats stuff, such as graph drawing.

use puruspe::{erf, gamma};
use std::f64::consts::PI;

const FRAC_1_SQRT_2PI: f64 = 0.3989422804014327; // (2.0f64 * PI).sqrt().recip()
const LN_SQRT_2PI: f64 = 0.9189385332046727; // (2.0 * PI).sqrt().ln()
const LN_SQRT_PID2: f64 = 0.22579135264472733; // (0.5 * PI).sqrt().ln()

pub fn dnorm(x: f64, mean: f64, sd: f64, log: bool) -> f64 {
    (sd * FRAC_1_SQRT_2PI).recip() * (-0.5 * ((x - mean) / sd).powi(2)).exp()
}

pub fn pnorm(q: f64, mean: f64, sd: f64, lower_tail: bool, log: bool) -> f64 {
    0.5 * (1.0 + erf((q - mean) / (sd * 2.0f64.sqrt())))
}

pub fn qnorm(p: f64, mean: f64, sd: f64, lower_tail: bool, log: bool) -> f64 {
    todo!()
}

/// Compute the probability density function for the t distribution at `x`. `df` is the degrees of
/// freedom and for hypothesis testing is the size of the sample - 1. `ncp` isn't implemented yet.
pub fn dt(x: f64, df: f64, ncp: Option<f64>, log: bool) -> f64 {
    let n = df;
    if x.is_nan() || !(n > 0.0) {
        return f64::NAN;
    }
    if ncp.is_some() {
        todo!()
    }
    if !x.is_finite() {
        return if log { f64::NEG_INFINITY } else { 0.0 };
    }
    if !n.is_finite() {
        return dnorm(x, 0.0, 1.0, log);
    }
    let t = -bd0(n * 0.5, (n + 1.0) * 0.5) + stirlerr((n + 1.0) * 0.5) - stirlerr(n * 0.5);
    let x2n = (x * x) / n;
    let lrg_x2n = x2n > f64::EPSILON.recip();
    let (l_x2n, u) = if lrg_x2n {
        let ax = x.abs();
        let l_x2n = ax.ln() - n.ln() * 0.5;
        (l_x2n, n * l_x2n)
    } else if x2n > 0.2 {
        let l_x2n = (1.0 + x2n).ln() * 0.5;
        (l_x2n, n * l_x2n)
    } else {
        let l_x2n = x2n.ln_1p() * 0.5;
        (l_x2n, -bd0(n * 0.5, (n + x * x) * 0.5) + x * x * 0.5)
    };
    if log {
        t - u - (LN_SQRT_2PI + l_x2n)
    } else {
        let i_sqrt = if lrg_x2n {
            n.sqrt() / x.abs()
        } else {
            (-l_x2n).exp()
        };
        (t - u).exp() * FRAC_1_SQRT_2PI * i_sqrt
    }
}

/// Compute the culmulative distribution function (CDF) for the t distribution at `x`. `df` is the
/// degrees of freedom and for hypothesis testing is the size of the sample - 1. `ncp`,
/// `lower_tail` and `log_p` aren't implemented yet.
pub fn pt(x: f64, df: f64, ncp: Option<f64>, lower_tail: bool, log_p: bool) -> f64 {
    assert!(df > 0.0, "degrees of freedom must be > 0");
    if ncp.is_some() || lower_tail || log_p {
        todo!()
    }
    todo!()
}

/// Compute the quantile function (inverse CDF function) for the t distribution at `x`.
pub fn qt(x: &[f64], df: f64, ncp: Option<f64>, lower_tail: bool, log_p: bool) -> f64 {
    assert!(df > 0.0, "degrees of freedom must be > 0");
    todo!()
}

fn bd0(x: f64, np: f64) -> f64 {
    if !x.is_finite() || !np.is_finite() || np == 0.0 {
        return f64::NAN;
    }

    if (x - np).abs() < 0.1 * (x + np) {
        let v = (x - np) / (x + np);
        let mut s = (x - np) * v;
        if s.abs() < f64::MIN_POSITIVE {
            return s;
        }
        let mut ej = 2.0 * x * v;
        let v2 = v * v;
        // prevents infinite loop, which is possible for certain input values.
        for j in 1..1000 {
            ej *= v2;
            let s1 = s + ej / ((j << 1) as f64 + 1.0);
            if s1 == s {
                return s1;
            }
            s = s1;
        }
    }

    // Either the loop didn't terminate or | x - np | is large enough for this to be accurate.
    x * (x / np).ln() + np - x
}

/// stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
///             = log Gamma(n+1) - 1/2 /// [log(2*pi) + log(n)] - n*[log(n) - 1]
///             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
///
/// see also lgammacor() in ./lgammacor.c  which computes almost the same!
fn stirlerr(n: f64) -> f64 {
    const S0: f64 = 0.083333333333333333333; // 1/12
    const S1: f64 = 0.00277777777777777777778; // 1/360
    const S2: f64 = 0.00079365079365079365079365; // 1/1260
    const S3: f64 = 0.000595238095238095238095238; // 1/1680
    const S4: f64 = 0.0008417508417508417508417508; // 1/1188

    // error for 0, 0.5, 1.0, 1.5, ..., 14.5, 15.0.
    const sferr_halves: &[f64] = &[
        0.0,                           /* n=0 - wrong, place holder only */
        0.1534264097200273452913848,   /* 0.5 */
        0.0810614667953272582196702,   /* 1.0 */
        0.0548141210519176538961390,   /* 1.5 */
        0.0413406959554092940938221,   /* 2.0 */
        0.03316287351993628748511048,  /* 2.5 */
        0.02767792568499833914878929,  /* 3.0 */
        0.02374616365629749597132920,  /* 3.5 */
        0.02079067210376509311152277,  /* 4.0 */
        0.01848845053267318523077934,  /* 4.5 */
        0.01664469118982119216319487,  /* 5.0 */
        0.01513497322191737887351255,  /* 5.5 */
        0.01387612882307074799874573,  /* 6.0 */
        0.01281046524292022692424986,  /* 6.5 */
        0.01189670994589177009505572,  /* 7.0 */
        0.01110455975820691732662991,  /* 7.5 */
        0.010411265261972096497478567, /* 8.0 */
        0.009799416126158803298389475, /* 8.5 */
        0.009255462182712732917728637, /* 9.0 */
        0.008768700134139385462952823, /* 9.5 */
        0.008330563433362871256469318, /* 10.0 */
        0.007934114564314020547248100, /* 10.5 */
        0.007573675487951840794972024, /* 11.0 */
        0.007244554301320383179543912, /* 11.5 */
        0.006942840107209529865664152, /* 12.0 */
        0.006665247032707682442354394, /* 12.5 */
        0.006408994188004207068439631, /* 13.0 */
        0.006171712263039457647532867, /* 13.5 */
        0.005951370112758847735624416, /* 14.0 */
        0.005746216513010115682023589, /* 14.5 */
        0.005554733551962801371038690, /* 15.0 */
    ];

    if n <= 15.0 {
        let nn = n + n;
        if nn == nn.trunc() {
            return sferr_halves[nn.trunc() as usize];
        }
        return lgammafn(n + 1.0, None) - (n + 0.5) * n.ln() + n - LN_SQRT_2PI;
    }

    let nn = (n * n).recip();
    if n > 500.0 {
        (S0 - S1 * nn) / n
    } else if n > 80.0 {
        (S0 - (S1 - S2 * nn) * nn) / n
    } else if n > 35.0 {
        (S0 - (S1 - (S2 - S3 * nn) * nn) * nn) / n
    } else {
        // 15 < n <= 35
        (S0 - (S1 - (S2 - (S3 - S4 * nn) * nn) * nn) * nn) / n
    }
}

/// First return is the result of ln(gamma(x)). Second is the sign of the gamma function: true
/// means positive and false means negative.
fn lgammafn(x: f64, sgn: Option<&mut bool>) -> f64 {
    const XMAX: f64 = 2.5327372760800758e+305; // f64::MAX / f64::MAX.ln()
    const DXREL: f64 = 1.4901161193847656e-08; // f64::EPSILON.sqrt()

    if x.is_nan() {
        return x;
    }
    if let Some(sgn) = sgn {
        *sgn = !(x < 0.0 && (-x).floor() % 2.0 == 0.0);
    }
    if x <= 0.0 && x == x.trunc() {
        return f64::INFINITY;
    }
    let y = x.abs();

    if y < 1e-306 {
        return -y.ln();
    }
    if y <= 10.0 {
        return gamma(x).abs().ln();
    }
    if y > XMAX {
        return f64::INFINITY;
    }
    if x > 0.0 {
        return if x > 1e17 {
            x * (x.ln() - 1.0)
        } else if x > 4934720.0 {
            LN_SQRT_2PI + (x - 0.5) * x.ln() - x
        } else {
            LN_SQRT_2PI + (x - 0.5) * x.ln() - x + lgammacor(x)
        };
    }
    let sinpiy = sinpi(y).abs();
    let ans = LN_SQRT_PID2;
    if ((x - (x - 0.5).trunc()) * ans / x) < DXREL {
        //eprintln!("warning lgammafn precision low as x near 0.");
    }
    ans
}

fn lgammacor(x: f64) -> f64 {
    const ALG_MCS: &[f64] = &[
        0.1666389480451863247205729650822e+0,
        -0.1384948176067563840732986059135e-4,
        0.9810825646924729426157171547487e-8,
        -0.1809129475572494194263306266719e-10,
        0.6221098041892605227126015543416e-13,
        -0.3399615005417721944303330599666e-15,
        0.2683181998482698748957538846666e-17,
        -0.2868042435334643284144622399999e-19,
        0.3962837061046434803679306666666e-21,
        -0.6831888753985766870111999999999e-23,
        0.1429227355942498147573333333333e-24,
        -0.3547598158101070547199999999999e-26,
        0.1025680058010470912000000000000e-27,
        -0.3401102254316748799999999999999e-29,
        0.1276642195630062933333333333333e-30,
    ];

    const NALGM: usize = 5;
    const XBIG: f64 = 94906265.62425156; // 2.0f64.powf(26.5)
    const XMAX: f64 = f64::MAX / 48.0;

    if x < 10.0 {
        f64::NAN
    } else if x >= XMAX {
        f64::NAN // would underflow
    } else if x < XBIG {
        let tmp = 10.0 / x;
        chebyshev_eval(tmp * tmp * 2.0 - 1.0, &ALG_MCS[..NALGM]) / x
    } else {
        (x * 12.0).recip()
    }
}

fn chebyshev_eval(x: f64, a: &[f64]) -> f64 {
    // sanity checks.
    if a.len() < 1 || a.len() > 1000 {
        return f64::NAN;
    }
    // I assume that the approximation is not accurate away from the origin.
    if x < -1.0 || x > 1.1 {
        return f64::NAN;
    }
    let xx = x + x;
    let mut b2 = 0.0;
    let mut b1 = 0.0;
    let mut b0 = 0.0;
    for coeff in a.into_iter().rev() {
        b2 = b1;
        b1 = b0;
        b0 = xx * b1 - b2 + coeff;
    }
    (b0 - b2) * 0.5
}

/// like x.sin(), except exact for multiples of 0.5.
fn sinpi(x: f64) -> f64 {
    let x = x % 2.0; // sin(pi(x + 2k)) == sin(pi x)  for all integer k
                     // map (-2,2) --> (-1,1] :
    let x = if x <= -1.0 {
        x + 2.0
    } else if x > 1.0 {
        x - 2.0
    } else {
        x
    };
    if x == 0.0 || x == 1.0 {
        0.0
    } else if x == 0.5 {
        1.0
    } else if x == -0.5 {
        -1.0
    } else {
        (PI * x).sin()
    }
}
