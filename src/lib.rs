//! This module provides functions from R's `stats` module, but rewritten for rust.
//!
//! I'm sticking to the `stats` module's naming for functions, even though they aren't at all
//! clear. Maybe this could change in the future.
//!
//! I'm also not implementing the non-stats stuff, such as graph drawing.

mod puruspe;

use puruspe::{betai, erf, gamma, inverf, ln_gamma};
use std::f64::consts::PI;

const FRAC_1_SQRT_2PI: f64 = 0.3989422804014327; // (2.0f64 * PI).sqrt().recip()
const LN_SQRT_2PI: f64 = 0.9189385332046727; // (2.0 * PI).sqrt().ln()
const LN_SQRT_PID2: f64 = 0.22579135264472733; // (0.5 * PI).sqrt().ln()
const LN_2PI: f64 = 1.8378770664093453; // (2.0 * PI).ln()

macro_rules! r_d__0 {
    ($log_p:expr) => {
        if $log_p {
            f64::NEG_INFINITY
        } else {
            0.0
        }
    };
}

macro_rules! r_d__1 {
    ($log_p:expr) => {
        if $log_p {
            0.0
        } else {
            1.0
        }
    };
}

macro_rules! r_dt_0 {
    ($lower_tail:expr, $log_p:expr) => {
        if $lower_tail {
            r_d__0!($log_p)
        } else {
            r_d__1!($log_p)
        }
    };
}

macro_rules! r_dt_1 {
    ($lower_tail:expr, $log_p:expr) => {
        if $lower_tail {
            r_d__1!($log_p)
        } else {
            r_d__0!($log_p)
        }
    };
}

macro_rules! r_d_cval {
    ($p:expr, $lower_tail:expr) => {
        if $lower_tail {
            0.5 - $p + 0.5
        } else {
            $p
        }
    };
}

macro_rules! r_d_exp {
    ($x:expr, $log:expr) => {
        if $log {
            $x
        } else {
            $x.exp()
        }
    };
}

macro_rules! r_d_fexp {
    ($f:expr, $x:expr, $log:expr) => {
        if $log {
            -0.5 * $f.ln() + $x
        } else {
            $x.ln() / $f.sqrt()
        }
    };
}

/// Evaluate the probability density function of the normal distribution at x.
pub fn dnorm(x: f64, mean: f64, sd: f64, log: bool) -> f64 {
    let out = (sd * FRAC_1_SQRT_2PI).recip() * (-0.5 * ((x - mean) / sd).powi(2)).exp();
    if log {
        out.ln()
    } else {
        out
    }
}

/// Evaluate the culmulative distribution function of the normal distribution at x.
pub fn pnorm(q: f64, mean: f64, sd: f64, lower_tail: bool, log: bool) -> f64 {
    let mut out = 0.5 * (1.0 + erf((q - mean) / (sd * 2.0f64.sqrt())));
    if !lower_tail {
        out = 1.0 - out;
    }
    if log {
        out.ln()
    } else {
        out
    }
}

/// Evaluate the quantile function of the normal distribution at probability p.
pub fn qnorm(mut p: f64, mean: f64, sd: f64, lower_tail: bool, log: bool) -> f64 {
    if !lower_tail {
        p = 1.0 - p;
    }
    let out = mean + sd * 2.0f64.sqrt() * inverf(2.0 * p - 1.0);
    if log {
        out.ln()
    } else {
        out
    }
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
/// degrees of freedom and for hypothesis testing is the size of the sample - 1.
pub fn pt(x: f64, df: f64, lower_tail: bool, log_p: bool) -> f64 {
    let n = df;
    if n <= 0.0 || n.is_nan() || x.is_nan() {
        return f64::NAN;
    }
    if x == f64::NEG_INFINITY {
        return r_dt_0!(lower_tail, log_p);
    } else if x == f64::INFINITY {
        return r_dt_1!(lower_tail, log_p);
    }
    if n == f64::INFINITY {
        return pnorm(x, 0.0, 1.0, lower_tail, log_p);
    }
    let nx = 1.0 + (x / n) * x;

    let val = if nx > 1e100 {
        let lval = -0.5 * n * (2.0 * x.abs().ln() - n.ln()) - lbeta(0.5 * n, 0.5) - (0.5 * n).ln();
        if log_p {
            lval
        } else {
            lval.exp()
        }
    } else {
        if n > x * x {
            pbeta(x * x / (n + x * x), 0.5, n * 0.5, false, log_p)
        } else {
            pbeta(nx.recip(), n * 0.5, 0.5, true, log_p)
        }
    };
    let lower_tail = if x <= 0.0 { !lower_tail } else { lower_tail };
    if log_p {
        if lower_tail {
            (-0.5f64 * val.exp()).ln_1p()
        } else {
            val - 2.0f64.ln()
        }
    } else {
        r_d_cval!(val * 0.5, lower_tail)
    }
}

// /// Compute the quantile function (inverse CDF function) for the t distribution at `x`.
// pub fn qt(x: f64, df: f64, ncp: Option<f64>, lower_tail: bool, log_p: bool) -> f64 {
//    if df <= 0.0 {
//        return f64::NAN;
//    }
//    todo!()
//}

//dbeta todo
/// Compute the probability density functino of the beta distribution at x.
pub fn pbeta(x: f64, a: f64, b: f64, lower_tail: bool, log_p: bool) -> f64 {
    if x.is_nan() || a.is_nan() || b.is_nan() {
        x + a + b
    } else if a < 0.0 || b < 0.0 {
        f64::NAN
    } else if x <= 0.0 {
        r_dt_0!(lower_tail, log_p)
    } else if x >= 1.0 {
        r_dt_1!(lower_tail, log_p)
    } else if a == 0.0 || b == 0.0 || !a.is_finite() || !b.is_finite() {
        if a == 0.0 && b == 0.0 {
            if log_p {
                -2.0f64.ln()
            } else {
                0.5
            }
        } else if a == 0.0 || a / b == 0.0 {
            r_dt_1!(lower_tail, log_p)
        } else if b == 0.0 || b / a == 0.0 {
            r_dt_0!(lower_tail, log_p)
        } else if x < 0.5 {
            r_dt_0!(lower_tail, log_p)
        } else {
            r_dt_1!(lower_tail, log_p)
        }
    } else {
        let mut w = betai(a, b, x);
        if !lower_tail {
            w = 1.0 - w;
        }
        if log_p {
            w = w.ln()
        }
        w
    }
}
//qbeta todo

/// Calculate the probability density function for the F distribution.
pub fn df(x: f64, m: f64, n: f64, log: bool) -> f64 {
    if x.is_nan() || m.is_nan() || n.is_nan() {
        x + m + n
    } else if m <= 0.0 || n <= 0.0 {
        f64::NAN
    } else if x < 0.0 {
        r_d__0!(log)
    } else if x == 0.0 {
        if m > 2.0 {
            r_d__0!(log)
        } else if m == 2.0 {
            r_d__1!(log)
        } else {
            f64::INFINITY
        }
    } else if m == f64::INFINITY && n == f64::INFINITY {
        if x == 1.0 {
            f64::INFINITY
        } else {
            r_d__0!(log)
        }
    } else if n == f64::INFINITY {
        dgamma(x, m / 2.0, 2.0 / m, log)
    } else if m > 1e14 {
        let dens = dgamma(1.0 / x, n / 2.0, 2.0 / n, log);
        if log {
            dens - 2.0 * x.ln()
        } else {
            dens / (x * x)
        }
    } else {
        let mut f = 1.0 / (n + x * m);
        let q = n * f;
        let p = x * m * f;
        let dens = if m >= 2.0 {
            f = m * q * 0.5;
            dbinom_raw((m - 2.0) * 0.5, (m + n - 2.0) * 0.5, p, q, log)
        } else {
            f = m * m * q / (2.0 * p * (m + n));
            dbinom_raw(m * 0.5, (m + n) * 0.5, p, q, log)
        };
        if log {
            f.ln() + dens
        } else {
            f * dens
        }
    }
}

// /// Calculates the culmulative distribution function for the F distribution at q.
// pub fn pf(q: f64, df1: f64, df2: f64, ncp: Option<f64>, lower_tail: bool, log_p: bool) -> f64 {
//    if ncp.is_some() {
//        todo!();
//    }
//    todo!()
//}
// qf

// dchisq
// pchisq
// qchisq

/// Calculates the probability density function for the gamma distribution at x.
pub fn dgamma(x: f64, shape: f64, scale: f64, log: bool) -> f64 {
    if x.is_nan() || shape.is_nan() || scale.is_nan() {
        x + shape + scale
    } else if shape < 0.0 || scale < 0.0 {
        f64::NAN
    } else if x < 0.0 {
        r_d__0!(log)
    } else if shape == 0.0 {
        if x == 0.0 {
            f64::INFINITY
        } else {
            r_d__0!(log)
        }
    } else if x == 0.0 {
        if shape < 1.0 {
            f64::INFINITY
        } else if shape > 1.0 {
            r_d__0!(log)
        } else if log {
            -scale.ln()
        } else {
            scale.recip()
        }
    } else if shape < 1.0 {
        let pr = dpois_raw(shape, x / scale, log);
        if log {
            pr + if (shape / x).is_finite() {
                (shape / x).ln()
            } else {
                // used only if (shape / x) overflows (to +inf).
                shape.ln() - x.ln()
            }
        } else {
            pr * shape / x
        }
    } else {
        let pr = dpois_raw(shape - 1.0, x / scale, log);
        if log {
            pr - scale.ln()
        } else {
            pr / scale
        }
    }
}

// pgamma
// qgamma

/// Evaluates the "deviance part" from Catherine Loader's algorithm.
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
    const SFERR_HALVES: &[f64] = &[
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
            return SFERR_HALVES[nn.trunc() as usize];
        }
        return lgammafn(n + 1.0) - (n + 0.5) * n.ln() + n - LN_SQRT_2PI;
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

/// Returns the result of ln(gamma(x)).
fn lgammafn(x: f64) -> f64 {
    const XMAX: f64 = 2.5327372760800758e+305; // f64::MAX / f64::MAX.ln()
    const DXREL: f64 = 1.4901161193847656e-08; // f64::EPSILON.sqrt()

    if x.is_nan() {
        return x;
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
    let ans = LN_SQRT_PID2 + (x - 0.5) * y.ln() - x - sinpiy.ln() - lgammacor(y);
    if ((x - (x - 0.5).trunc()) * ans / x) < DXREL {
        //eprintln!("warning lgammafn precision low as x near 0.");
    }
    ans
}

/// Compute the log gamma correction factor for x >= 10 so that
///
/// `log(gamma(x)) = .5*log(2*pi) + (x-.5)*log(x) -x + lgammacor(x)`
///
/// `[ lgammacor(x) is called Del(x) in other contexts (e.g. dcdflib)]
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

fn gammafn(x: f64) -> f64 {
    const GMACS: &[f64] = &[
        0.8571195590989331421920062399942e-2,
        0.4415381324841006757191315771652e-2,
        0.5685043681599363378632664588789e-1,
        -0.4219835396418560501012500186624e-2,
        0.1326808181212460220584006796352e-2,
        -0.1893024529798880432523947023886e-3,
        0.3606925327441245256578082217225e-4,
        -0.6056761904460864218485548290365e-5,
        0.1055829546302283344731823509093e-5,
        -0.1811967365542384048291855891166e-6,
        0.3117724964715322277790254593169e-7,
        -0.5354219639019687140874081024347e-8,
        0.9193275519859588946887786825940e-9,
        -0.1577941280288339761767423273953e-9,
        0.2707980622934954543266540433089e-10,
        -0.4646818653825730144081661058933e-11,
        0.7973350192007419656460767175359e-12,
        -0.1368078209830916025799499172309e-12,
        0.2347319486563800657233471771688e-13,
        -0.4027432614949066932766570534699e-14,
        0.6910051747372100912138336975257e-15,
        -0.1185584500221992907052387126192e-15,
        0.2034148542496373955201026051932e-16,
        -0.3490054341717405849274012949108e-17,
        0.5987993856485305567135051066026e-18,
        -0.1027378057872228074490069778431e-18,
        0.1762702816060529824942759660748e-19,
        -0.3024320653735306260958772112042e-20,
        0.5188914660218397839717833550506e-21,
        -0.8902770842456576692449251601066e-22,
        0.1527474068493342602274596891306e-22,
        -0.2620731256187362900257328332799e-23,
        0.4496464047830538670331046570666e-24,
        -0.7714712731336877911703901525333e-25,
        0.1323635453126044036486572714666e-25,
        -0.2270999412942928816702313813333e-26,
        0.3896418998003991449320816639999e-27,
        -0.6685198115125953327792127999999e-28,
        0.1146998663140024384347613866666e-28,
        -0.1967938586345134677295103999999e-29,
        0.3376448816585338090334890666666e-30,
        -0.5793070335782135784625493333333e-31,
    ];

    const NGAM: usize = 22;
    const XMIN: f64 = -170.5674972726612;
    const XMAX: f64 = 171.61447887182298;
    const XSML: f64 = 2.2474362225598545e-308;
    const DXREL: f64 = 1.490116119384765696e-8;
    if x.is_nan() {
        return x;
    }
    if x == 0.0 || x < 0.0 && x == x.round() {
        return f64::NAN;
    }

    let y = x.abs();
    if y <= 10.0 {
        let mut n = x as isize;
        if x < 0.0 {
            n -= 1;
        }
        let y = x - n as f64;
        n -= 1;
        let mut value = chebyshev_eval(y * 2.0 - 1.0, &GMACS[..NGAM]) + 0.9375;
        if n == 0 {
            return value;
        } else if n < 0 {
            if x < -0.5 && (x - (x - 0.5).floor() / x) < DXREL {
                //panic!("less than 1/2 precision as near a pole (negative integer)
            }
            if y < XSML {
                // overflow
                if x > 0.0 {
                    f64::INFINITY
                } else {
                    f64::NEG_INFINITY
                }
            } else {
                let n = -n as usize;
                for i in 0..n {
                    value /= x + i as f64;
                }
                value
            }
        } else {
            for i in 1..=n {
                value *= y + i as f64;
            }
            value
        }
    } else {
        if x > XMAX {
            return f64::INFINITY;
        } else if x < XMIN {
            return 0.0;
        }

        let mut value;
        if y <= 50.0 && y.floor() == y {
            value = 1.0;
            for i in 2..(y as usize) {
                value *= i as f64;
            }
        } else {
            let cor = if (2.0 * y).floor() == 2.0 * y {
                stirlerr(y)
            } else {
                lgammacor(y)
            };
            value = ((y - 0.5) * y.ln() - y + LN_SQRT_2PI + cor).exp();
        }

        if x > 0.0 {
            return value;
        }

        if (x - (x - 0.5).floor() / x).abs() < DXREL {
            //println!("less than 0.5 precision")
        }

        let sinpiy = sinpi(y);
        if sinpiy == 0.0 {
            f64::INFINITY
        } else {
            -PI / (y * sinpiy * value)
        }
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

/// Computes ln(beta(a, b))
fn lbeta(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return f64::NAN;
    }
    let p = a.min(b);
    let q = a.max(b);

    if p < 0.0 {
        // also checks q >= 0
        return f64::NAN;
    } else if p == f64::INFINITY {
        return f64::INFINITY;
    } else if q == f64::INFINITY {
        return f64::NEG_INFINITY;
    }

    if p >= 10.0 {
        let corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
        q.ln() * -0.5
            + LN_SQRT_2PI
            + corr
            + (p - 0.5) * (p / (p + q)).ln()
            + q * (-p / (p + q)).ln_1p()
    } else if q >= 10.0 {
        let corr = lgammacor(q) - lgammacor(p + q);
        lgammafn(p) + corr + p - p * (p + q).ln() + (q - 0.5) * (-p / (p + q)).ln_1p()
    } else if p < 1e-306 {
        ln_gamma(p) + (ln_gamma(q) - ln_gamma(p + q))
    } else {
        (gammafn(p) * (gammafn(q) / gammafn(p + q))).ln()
    }
}

fn dbinom_raw(x: f64, n: f64, p: f64, q: f64, log: bool) -> f64 {
    if p == 0.0 {
        if x == 0.0 {
            r_d__1!(log)
        } else {
            r_d__0!(log)
        }
    } else if q == 0.0 {
        if x == n {
            r_d__1!(log)
        } else {
            r_d__0!(log)
        }
    } else if x == 0.0 {
        if n == 0.0 {
            r_d__1!(log)
        } else {
            r_d_exp!(
                if p < 0.1 {
                    -bd0(n, n * q) - n * p
                } else {
                    n * q.ln()
                },
                log
            )
        }
    } else if x == n {
        r_d_exp!(
            if q < 0.1 {
                -bd0(n, n * p) - n * q
            } else {
                n * p.ln()
            },
            log
        )
    } else if x < 0.0 || x > n {
        r_d__0!(log)
    } else {
        let lc = stirlerr(n) - stirlerr(x) - stirlerr(n - x) - bd0(x, n * p) - bd0(n - x, n * q);
        let lf = LN_2PI + x.ln() + (-x / n).ln_1p();
        r_d_exp!(lc - 0.5 * lf, log)
    }
}

fn dpois_raw(x: f64, lambda: f64, log: bool) -> f64 {
    debug_assert!(x >= 0.0);
    debug_assert!(lambda >= 0.0);
    if lambda == 0.0 {
        if x == 0.0 {
            r_d__1!(log)
        } else {
            r_d__0!(log)
        }
    } else if lambda == f64::INFINITY {
        r_d__0!(log)
    } else if x < 0.0 {
        r_d__0!(log)
    } else if x <= lambda + f64::MIN_POSITIVE {
        r_d_exp!(-lambda, log)
    } else if lambda < x * f64::MIN_POSITIVE {
        if x == f64::INFINITY {
            r_d__0!(log)
        } else {
            r_d_exp!(-lambda + x * lambda.ln() - lgammafn(x + 1.0), log)
        }
    } else {
        r_d_fexp!(PI * 2.0 * x, -stirlerr(x) - bd0(x, lambda), log)
    }
}
