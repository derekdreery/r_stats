use crate::LN2;

#[inline]
#[allow(non_snake_case)]
pub fn r_d__0(log_p: bool) -> f64 {
    if log_p {
        f64::NEG_INFINITY
    } else {
        0.0
    }
}

#[inline]
#[allow(non_snake_case)]
pub fn r_d__1(log_p: bool) -> f64 {
    if log_p {
        0.0
    } else {
        1.0
    }
}

#[inline]
pub fn r_dt_0(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        r_d__0(log_p)
    } else {
        r_d__1(log_p)
    }
}

#[inline]
pub fn r_dt_1(lower_tail: bool, log_p: bool) -> f64 {
    if lower_tail {
        r_d__1(log_p)
    } else {
        r_d__0(log_p)
    }
}

#[inline]
pub fn r_d_cval(p: f64, lower_tail: bool) -> f64 {
    if lower_tail {
        0.5 - p + 0.5
    } else {
        p
    }
}

#[inline]
pub fn r_d_exp(x: f64, log: bool) -> f64 {
    if log {
        x
    } else {
        x.exp()
    }
}

#[inline]
pub fn r_d_fexp(f: f64, x: f64, log: bool) -> f64 {
    if log {
        -0.5 * f.ln() + x
    } else {
        x.ln() / f.sqrt()
    }
}

#[inline]
pub fn r_log1_exp(x: f64) -> f64 {
    if x > -LN2 {
        (-x.exp_m1()).ln()
    } else {
        (-x.exp()).ln_1p()
    }
}
