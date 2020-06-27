fn main() {
    cc::Build::new()
        .define("MATHLIB_STANDALONE", "1")
        .include("nmath")
        .files(&[
            "nmath/mlutils.c",
            "nmath/fmax2.c",
            "nmath/fmin2.c",
            "nmath/fprec.c",
            "nmath/fround.c",
            "nmath/ftrunc.c",
            "nmath/sign.c",
            "nmath/fsign.c",
            "nmath/imax2.c",
            "nmath/imin2.c",
            "nmath/chebyshev.c",
            "nmath/lgammacor.c",
            "nmath/gammalims.c",
            "nmath/stirlerr.c",
            "nmath/bd0.c",
        ])
        .compile("nmath")
}
