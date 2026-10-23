# Documentation

Start with the [project README](../README.md) for commands and output paths.

| Guide | Contents |
| --- | --- |
| [Data](data.md) | Tree layout, campaign splits, joins, units and missing values |
| [Photons](photons.md) | Selections, matching, detector regions and efficiency definitions |
| [BDT](bdt.md) | Feature sets, class IDs and training defaults |
| [Converter](converter.md) | Gamma/RICH fixes and remaining limitations |
| [Native branches](branches/native.md) | Historical Florian scan of `nanoaod.root` |
| [Raw branches](branches/raw.md) | Historical Florian scan of `nanoaod_raw_sdst.root` |
| [Classifier inputs](branches/classification.md) | Current porting checks for all 150 supplied inputs |

The Florian files predate the converter fixes; the Jongwon files use the
corrected Gamma/RICH schema. Updating the code does not update existing ROOT files.

Numeric scan tables cover the original Florian files:
[native](scans/native.csv) and [raw](scans/raw.csv).
