# Inputs and Outputs for Radiocarbon Calibration Pipeline

They live here.

## Tutorials
Here, we present a couple of usage scenarios.  
Before any of the steps below, you should:
1. Export the data to a CSV file
2. Make sure the column that has the C14 dates is called `C14Age`, and the column that has the error is called `C14SD`

### Scenario 1: New file to calibrate 
Let's say the file that we want to calibrate is called `example.csv`:

1. Open this repository and navigate to the `inputs` folder
2. Add the new CSV file, `example.csv`, for calibration:
   1. Top right corner
   2. `Add file`
   3. `Upload files`
   4. `Choose your files`
   5. Choose the file you want to calibrate
   6. `Commit changes`
3. Watch the status of the pipeline change from `In progress` to `Success` or `Failure`
4. In case of `Success`, check the `outputs` folder for a file called `example.radiocarbon.csv` and for new PDFs, those are your outputs! 

### Scenario 2: Run a file with another script and confidence value 
Let's say the file that we want to calibrate is called `sedimentary_charcoal.csv`:

1. Open this repository and navigate to the `inputs` folder
2. Create a new folder file inside `inputs`, inside the folder we will create a `config` file:
    1. Top right corner
    2. `Add file`
    3. `Create new file`
    4. Write "/sediment_analysis/config"  
       (The folder name can be whatever you want!)
    5. Copy the structure of the [default config file](./config) to the text field
    6. Change confidence to the desired value, and script to the name of a script that exists on the `scripts` folder, for example `sediments.r`
    7. On the top right, `Commit changes...`
    8. `Commit changes`  
3. Navigate to the new folder you just created
4. Add the new CSV file, `sedimentary_charcoal.csv`, for calibration (Top right corner -> `Add file` -> `Upload files` -> `Choose your files` -> Choose the file you want to calibrate -> `Commit changes`)
5. Watch the status of the pipeline change from `In progress` to `Success` or `Failure`
6. In case of `Success`, check the `outputs` folder, then `sediment_analysis` for a file called `sedimentary_charcoal.sediments.csv` and for new PDFs, those are your outputs!

---

## `config` files
You can modify the `config` files to adjust settings such as confidence intervals, time steps, filtering, and the script to be used for a file.
For example:
```file
step=20
confidence=0.90
column=Site
value=SMM
script=my_script.r
```
This configuration uses a 20-year step, a confidence interval of 90%, filter on all rows with the value `SMM` on the column named `Site`, and runs the script `my_script.r` on the file.

> [!WARNING]
> The order of values is important!  
> All values can be empty, except for `script`.

### Default and folder `config`s
Folder-specific `config` files override the default for all CSV files within that folder. If a file lacks a folder `config`, the default `config` is applied. Folder `config` files have the same structure and content as the default `config`.

By organizing your scripts and configurations this way, you can easily manage and customize the processing of different input files.

> [!WARNING]
> Whenever a folder `config` file is changed, the outputs for all files in that folder get recomputed!

---

## Folder Structure
The structure of the directories looks like this:
```
.
├── config                         # Default configuration for input files. Files without a nearby config use this.
├── inputs/                        # Location of input files.
│   ├── folder_1/                  # Files can be organized in folders.
│   │   ├── example1.csv           # This file uses the config next to it.
│   │   └── config
│   ├── folder_2/
│   │   └── example2.csv           # This file uses the default config.
│   └── example3.csv               # And this one too.
├── outputs/                       # Output files are stored here. Folder structure mirrors the inputs.
│   ├── folder_1/
│   │   ├── example1.my_script.csv # Processed file name includes the script used.
│   │   └── example1.pdf           # PDF files generated by the script are also stored here.
│   ├── folder_2/
│   │   └── example2.default_script.csv
│   └── example3.default_script.csv
└── scripts/                       # Location of the scripts.
    ├── default_script.r
    └── my_script.r
```
You can nest even more folders inside `inputs` if you like! It's all about organization. The structure is automatically mirrored on the outputs, so everything stays where you expect it to.

---

## Scripts 
For more information about how to create more scripts, check out [their information page](./scripts/README.md).
