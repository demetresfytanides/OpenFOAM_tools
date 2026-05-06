# mappedFlowRateNonUniformCorrected

A variant of `mappedFlowRateCorrected` that enforces a target water volumetric flow rate
by **multiplicatively rescaling** the normal velocity component instead of applying a
**uniform additive shift**.

## Why this BC exists

The original `mappedFlowRateCorrected` writes back

```
U[face] = Ut + (Un + alpha * dUn) * n,   dUn = (-Qtarget - QwRec)/Awet
```

`dUn` is a single scalar shared by every face on the patch, so it shifts the *mean* of the
normal component while preserving the absolute size of fluctuations. If `dUn` is large
compared to the rms of the donor's normal component, the BC distorts the donor's profile
shape and turbulence intensity ratios.

This variant writes back

```
U[face] = Ut + gammaEff(face) * Un * n,
gamma     = -Qtarget / QwRec   (clipped to [gammaMin, gammaMax])
gammaEff  = (1 - alpha) + alpha * gamma
```

so the donor's normal-component profile is rescaled by a single factor `gamma`. Mean and
fluctuations are scaled proportionally, ratios are preserved, and dry faces (`alpha = 0`)
receive the donor field untouched.

Tangential donor structure is passed through unchanged, exactly as in
`mappedFlowRateCorrected`.

## Files

```
mappedFlowRateNonUniformCorrectedFvPatchField.H
mappedFlowRateNonUniformCorrectedFvPatchField.C
mappedFlowRateNonUniformCorrectedFvPatchFields.H
mappedFlowRateNonUniformCorrectedFvPatchFields.C
Make/files
Make/options
```

The mapping infrastructure (`mappedFlowRateCorrectedPatchFieldBase.H/.C`) is reused from
the sibling library — `Make/options` adds `-I../mappedFlowRateCorrected` so the headers
are found.

## Build

```bash
cd $WM_PROJECT_USER_DIR/src/mappedFlowRateNonUniformCorrected
wclean
wmake libso
ls -l $FOAM_USER_LIBBIN/libmappedFlowRateNonUniformCorrected*
```

Then in `system/controlDict`:

```foam
libs
(
    "libmappedFlowRateNonUniformCorrected.so"
);
```

If you want both BCs available in the same case, also load `libmappedFlowRateCorrected.so`.

## Dictionary entry

```foam
inlet
{
    type                mappedFlowRateNonUniformCorrected;

    field               U;
    setAverage          no;
    average             uniform (0 0 0);
    interpolationScheme cell;

    Qtarget             0.089;       // target water volumetric flow rate
    alphaMin            1e-3;        // wet-cell threshold
    gammaMin            0.5;         // lower clip on rescale factor
    gammaMax            2.0;         // upper clip on rescale factor

    sampleMode          nearestCell;
    sampleRegion        region0;
    samplePatch         none;
    offsetMode          uniform;
    offset              (0 2 0);

    value               uniform (0 0 0);
}
```

## Tuning notes

- `gammaMin` / `gammaMax` are safety clips. In a healthy run `gamma` should sit close to
  `1.0`; values that hug the clips persistently mean the donor's water Q is far from
  `Qtarget` and you should re-check the donor.
- `alphaMin` excludes nearly-dry faces from `QwRec` and `Awet`. Keep it small enough that
  legitimately wet faces are not skipped (`1e-3` to `1e-2` is typical).
- Run with `debug 1` for this BC class in `controlDict`'s `DebugSwitches` to print
  `gamma`, `QwRec`, and post-correction `QwAfter` each step.

## Validation checklist

After running, verify:

```text
Qw_post  = sum( alpha * (U . n) * Sf )   ≈ -Qtarget
gamma    ∈ [gammaMin, gammaMax]          and ideally near 1.0
```
