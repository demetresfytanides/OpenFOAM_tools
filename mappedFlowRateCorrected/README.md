# README — mappedFlowRateCorrected on Improv

## What this branch adds

This branch adds a custom OpenFOAM boundary condition named `mappedFlowRateCorrected`.

The goal is to use a donor plane to map inlet velocity structure while enforcing a target **water** volumetric flow rate:

- the donor velocity field is mapped to the inlet using OpenFOAM's mapped-field machinery
- local `alpha.water` at the inlet is used to identify the water fraction on each face
- only the **normal** component of velocity is corrected
- the **tangential** components are preserved from the donor field
- the target volumetric water flow rate `Qtarget` is read from `0/U`

In words, the boundary condition does this each time it updates:

1. maps the donor velocity field onto the inlet
2. computes the current inlet water flux
3. computes the wet area using local `alpha.water`
4. computes a scalar normal correction `dUn`
5. applies that correction only to the normal velocity component, weighted by `alpha.water`

This preserves the donor turbulence structure as much as possible while still enforcing the desired water discharge.

## Files added or modified

The custom library lives under:

```text
$WM_PROJECT_USER_DIR/src/mappedFlowRateCorrected
```

Important files:

```text
mappedFlowRateCorrectedFvPatchField.H
mappedFlowRateCorrectedFvPatchField.C
mappedFlowRateCorrectedFvPatchFields.H
mappedFlowRateCorrectedFvPatchFields.C
mappedFlowRateCorrectedPatchFieldBase.H
mappedFlowRateCorrectedPatchFieldBase.C
Make/files
Make/options
```

## Expected runtime dictionary entries

The inlet entry in `0/U` should look like this:

```foam
inlet
{
    type                mappedFlowRateCorrected;

    field               U;
    setAverage          no;
    average             uniform (0 0 0);
    interpolationScheme cell;

    Qtarget             0.089;
    alphaMin            1e-8;

    sampleMode          nearestCell;
    sampleRegion        region0;
    samplePatch         none;
    offsetMode          uniform;
    offset              (0 2 0);

    value               uniform (0 0 0);
}
```

If restarting from a nonzero time, the same inlet entry must also exist in the restart directory, for example `1350/U`.

## Git workflow

### 1. Create a feature branch locally

```bash
git checkout -b feature/mappedFlowRateCorrected
```

### 2. Add the new source files and case changes

```bash
git add   $WM_PROJECT_USER_DIR/src/mappedFlowRateCorrected   system/controlDict   0/U
```

### 3. Commit

```bash
git commit -m "Add mappedFlowRateCorrected inlet BC for mapped inflow with alpha-weighted normal flux correction"
```

### 4. Push the new branch

```bash
git push -u origin feature/mappedFlowRateCorrected
```

## Accessing Improv

Improv access is through the CELS login nodes using a jump-host workflow. Do not assume direct SSH to Improv is allowed.

## Recommended directory layout on Improv

After logging in, work in your project filesystem, not in your home directory.

Example:

```bash
cd /path/to/your/project
mkdir -p openfoam-dev
cd openfoam-dev
git clone <your-repo-url>
cd <your-repo-name>
git checkout feature/mappedFlowRateCorrected
```

## Sourcing the environment

### Option A — site-provided OpenFOAM installation

If Improv provides an OpenFOAM build, use that installation's `etc/bashrc`.

```bash
source /path/to/OpenFOAM-v2506/etc/bashrc
```

### Option B — your own OpenFOAM installation

```bash
source /path/to/your/OpenFOAM-v2506/etc/bashrc
```

### Verify the environment

```bash
echo $WM_PROJECT_VERSION
echo $WM_PROJECT_DIR
echo $FOAM_USER_LIBBIN
which wmake
```

### Helpful `.bashrc` block

```bash
if [ -f /etc/bashrc ]; then
    . /etc/bashrc
elif [ -f /etc/bash.bashrc ]; then
    . /etc/bash.bashrc
fi

source /path/to/your/OpenFOAM-v2506/etc/bashrc
```

Then reload:

```bash
source ~/.bashrc
```

## Building the custom BC on Improv

```bash
cd $WM_PROJECT_USER_DIR/src/mappedFlowRateCorrected
wclean
wmake libso
ls -l $FOAM_USER_LIBBIN/libmappedFlowRateCorrected*
```

## Load the library in the case

In `system/controlDict`, add:

```foam
libs
(
    "libmappedFlowRateCorrected.so"
);
```

## Running the case

Before running, check that:
- `0/U` or the restart-time `U` file contains the `mappedFlowRateCorrected` inlet entry
- `system/controlDict` loads the library
- the mapped sampling parameters are correct
- `Qtarget` is set correctly

Then run normally, for example:

```bash
interFoam > log.interFoam 2>&1
```

## What to verify in the results

Visualize the inlet and donor plane for:
- `Uy`
- `Ux`
- `Uz`

Check that the inlet water flux is close to `Qtarget`, where:

```text
Qw = Σ alpha.water * (U · n) * A
```

## Suggested command sequence on Improv

```bash
cd /path/to/project
mkdir -p openfoam-dev
cd openfoam-dev

git clone <your-repo-url>
cd <your-repo-name>
git checkout feature/mappedFlowRateCorrected

source /path/to/OpenFOAM-v2506/etc/bashrc

echo $WM_PROJECT_VERSION
which wmake

cd $WM_PROJECT_USER_DIR/src/mappedFlowRateCorrected
wclean
wmake libso

cd /path/to/case
grep -n "mappedFlowRateCorrected" 0/U system/controlDict
interFoam > log.interFoam 2>&1
```

## Final note

This BC is intended specifically for velocity (`U`) and was compiled vector-only on purpose. It is not meant to be a generic scalar/tensor mapped condition.
