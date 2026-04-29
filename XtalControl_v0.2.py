import pandas as pd
import numpy as np
from natsort import natsort_keygen

metadata = {
    "Author" : "Milo R. Cooper",
    "Version" : "v0.2",
    "Program" : "XtalControl",
    }

# ----------- Pulling in data / trimming -----------------------

xcePath = 'python_scripts/XtalControl/Raw_data/3C_xceexport_marchsometime.csv'
xceImport = pd.read_csv(xcePath)

codePath = 'python_scripts/XtalControl/control-plate-map.csv'
codeImport = pd.read_csv(codePath)                 # This is used to trim the rows in xceImport

compoundStockConcentration = 500 # mM !! Check this !!

possibleControls = ["DMSO","DMSO30", "30", "DMSO25", "25", "DMSO20", "20", "DMSO15", "15", "DMSO10", "10", "DMSO5", "5"] # may need updating as people find other stupid ways to refer to DMSO
controlCodes = pd.DataFrame({"Code" : possibleControls})
codeImport = pd.concat([codeImport, controlCodes], ignore_index = True)

xceImport = xceImport.copy()
xceImport["KeepCode"] = (xceImport.CompoundCode.isin(codeImport.Code)) # Boolean column for whether a code is present in the control platemap + possibleControls
xceImport = xceImport[xceImport.KeepCode] # drops rows where KeepCode = False
xceImport.drop(columns="KeepCode") # drops KeepCode (no longer needed)

# ------------ FAILURE PARAMETERS -------------------

mountCutoffFailure = 49
mountCutoffPartial = 70

diffractCutoffFailure = 0.1                # this just makes the parameter editable, not expected to change: (diff attempt - diff success) > 0.1 = FAIL

resolutionCutoffFailure = 3.0
resolutionCutoffPartial = 2.2

dimpleRunCutoffFailure = 34                # (Diffract success - Dimple success) < cutoff = Pass 
dimpleRunCutoffPartial = 0

rCrystCutoff = 0.40                        # Cutoff boundary
rFreeCutoff = 0.35

modelBuildingCutoffFailure = 49            # i.e. diffract attempt % - Dimple success % > 49 = FAIL (100 - 33 = 66, 66 > 49 = Fail)
modelBuildingCutoffPartial = 25            # e.g. 100 - 66 = 33, 33 > 25 = PARTIAL (fail rate of 1 in 3)

# -------------- END PARAMETERS -----------------------

splitCols = (
    xceImport['MountingResult']
    .str.split(':', n=2, expand=True)
    .set_axis(['MountResult', 'DropComment', 'CrystalComment'], axis=1)
    .apply(lambda x: x.str.strip())
)                          
xceImport = pd.concat([xceImport, splitCols], axis=1)                                       # Creates 3 new columns and trims whitespace from new columns - no lambda throws type error
xceImport.drop(columns=['MountingResult'], inplace=True)                                    # deletes old column to clean up

xceImport = xceImport.copy()                                                                # fixes frag error, which shouldn't be happing anyway as I'm using concat.
identCodes = (xceImport["CompoundCode"] + "/" + xceImport["SolventFraction"].astype(str))   
xceImport["SoakCodes"] = identCodes                                                         # Create a single-column identifying code (e.g. "Nitromethane/10") that we can divide data by from now on - maybe add soaktime here?

unique_soaks = xceImport["SoakCodes"].unique()
soaksCollated = pd.DataFrame({"SoakCodes": unique_soaks})                                   # soaksCollated will be the main dataframe, averaged across unique experiments.

counted_soaks = xceImport.value_counts("SoakCodes")                                         
soaksCollated = pd.merge(soaksCollated, counted_soaks, how = 'left', on = 'SoakCodes', validate = "1:1")
soaksCollated = soaksCollated.rename(columns = {"count": "SoaksAttempted"})
soaksCollated = soaksCollated.sort_values('SoakCodes')                                      # counts incidence of each code = number of soaks attempted

def reportIncidence(targetColumn):
    repincidence = (                                    # Reads through target column and counts the values as %, grouping by SoakCode + unstacks the produced table into multiple columns
        xceImport.groupby("SoakCodes")[targetColumn]                                           
        .value_counts(normalize=True, dropna=False)
        .mul(100)
        .round(1)
        .unstack(fill_value=0)
    )
    return repincidence

mountSuccess = reportIncidence("MountResult")
mountSuccess = mountSuccess.rename(columns = {"FAIL": "Mount-FAIL %", "OK": "Mount-OK %"})
mountSuccess = mountSuccess.sort_values('SoakCodes')
soaksCollated = pd.merge(soaksCollated, mountSuccess, how = 'left', on = 'SoakCodes', validate = "1:1")
soaksCollated.drop(columns=[''], inplace=True, errors = "ignore")  # drops blank column, error ignored in case there's no blank cells

dropComment = reportIncidence("DropComment")
dropComment = dropComment.rename(columns = {"Precipitated": "Drop-Precipitated %", "Crystalline": "Drop-Crystalline %", "Bad Dispense": "Drop-Bad Dispense %"})
soaksCollated = pd.merge(soaksCollated, dropComment, how = 'left', on = 'SoakCodes', validate = "1:1")
soaksCollated.drop(columns=[''], inplace=True, errors = "ignore")                                                            

crystalComment = reportIncidence("CrystalComment")
crystalComment = crystalComment.rename(columns = {"Cracked": "Crystal-Cracked %", "Dissolved": "Crystal-Dissolved %", "Jelly": "Crystal-Jelly %"})
soaksCollated = pd.merge(soaksCollated, crystalComment, how = 'left', on = 'SoakCodes', validate = "1:1")
soaksCollated.drop(columns=[''], inplace=True, errors = "ignore")

all_soakcodes = xceImport.value_counts('SoakCodes')                                     # counts row per soakcode in full dataframe

soaksCollated["Diffraction Attempt %"] = soaksCollated['SoakCodes'].map(                # uses .map because soaksCollated is integer indexed (change?)         
    xceImport[xceImport['CrystalName'].replace('', pd.NA).notna()]                      # Filters out rows with blank/single-space CrystalName
    .value_counts('SoakCodes')
    .reindex(all_soakcodes.index, fill_value=0)                                         # adds a 0 for any NaN value_counts where no soakcode repeat made it to a crystal name
    .div(all_soakcodes)
    .mul(100)
    .round(2)
)  # !! this is the same as 'Mount:OK %' because XCE doesn't log beamline data collection attempts, only successfully indexed diffraction datasets
   # !! so ISPyB submission is the last thing logged before a successful data collection is registered. ISPyB sample csv doesn't log either.

soaksCollated["Diffraction Success %"] = soaksCollated['SoakCodes'].map(  # this should be a function really         
    xceImport[xceImport['DataCollectionOutcome'].replace('', pd.NA).notna()]                      
    .value_counts('SoakCodes')
    .reindex(all_soakcodes.index, fill_value=0)                                         
    .div(all_soakcodes)
    .mul(100)
    .round(2)
)     

codeResolution = xceImport.groupby("SoakCodes")["DataProcessingResolutionHigh"].mean(skipna=True, numeric_only=True).round(2).fillna("NO DATA")
soaksCollated = (pd
    .merge(soaksCollated, codeResolution, how = 'left', on = 'SoakCodes', validate = "1:1")
    .rename(columns={"DataProcessingResolutionHigh":"Resolution (A)"})
)    

DIMPLESuccess = reportIncidence("DataProcessingDimpleSuccessful")
DIMPLESuccess = DIMPLESuccess.rename(columns = {True: 'DIMPLE Success %'}).sort_values('SoakCodes')
soaksCollated = pd.merge(soaksCollated, DIMPLESuccess, how = 'left', on = 'SoakCodes', validate = "1:1")
soaksCollated.drop(columns=np.nan, inplace=True, errors = "ignore")

codeRcryst = xceImport.groupby("SoakCodes")["DimpleRcryst"].mean(skipna=True, numeric_only=True).round(2).fillna("NO DATA")
soaksCollated = (pd
    .merge(soaksCollated, codeRcryst, how = 'left', on = 'SoakCodes', validate = "1:1")
    .rename(columns={"DimpleRcryst":"DIMPLE Rcryst"})
)

codeRfree = xceImport.groupby("SoakCodes")["DimpleRfree"].mean(skipna=True, numeric_only=True).round(2).fillna("NO DATA")  # All average data now pulled in, assess cutoffs?
soaksCollated = (pd
    .merge(soaksCollated, codeRfree, how = 'left', on = 'SoakCodes', validate = "1:1")
    .rename(columns={"DimpleRfree":"DIMPLE Rfree"})
)

def calcRPassRate(targetColumn, targetCutoff):
    numericDimpleRvalue = pd.to_numeric(xceImport[targetColumn], errors="coerce")  # forces numeric so blanks are NaN
    hasRData = numericDimpleRvalue.notna().groupby(xceImport["SoakCodes"]).any()  # checks which soakcodes have at least one Rcryst value (non-blank)
    passingCounts = xceImport[numericDimpleRvalue < targetCutoff].value_counts("SoakCodes")     
    totalCounts = xceImport.value_counts("SoakCodes")                                           
    rPassRate = passingCounts.reindex(totalCounts.index, fill_value=0).div(totalCounts).mul(100).round(2)
    rPassRate = rPassRate.astype(object)
    rPassRate[~hasRData] = "NO DATA"  # End result: NaN = "NO DATA"; no passes in the value_count = 0.00
    return rPassRate                                                      
    
rCrystRate = calcRPassRate("DimpleRcryst", rCrystCutoff)                                                
soaksCollated["Rcryst Pass %"] = soaksCollated["SoakCodes"].map(rCrystRate)  # These values don't account for blank values so if 1/3 crystals get a model and that model passes the cutoff, the rate is 33% - maybe change?

rFreeRate = calcRPassRate("DimpleRfree", rFreeCutoff)
soaksCollated["Rfree Pass %"] = soaksCollated["SoakCodes"].map(rFreeRate)

# --------------------- Start of parameter comparisons and grading ---------------------

mountChoices = ["FAIL", "PARTIAL"]
mountConditions = [
    soaksCollated["Mount-OK %"] <= mountCutoffFailure,
    soaksCollated["Mount-OK %"] <= mountCutoffPartial,
]
soaksCollated["Mount Grade"] = np.select(mountConditions, mountChoices, default = "PASS")  # First mountcondition assigned FAIL, second assigned PARTIAL, otherwise PASS

diffractionChoices = ["NO DATA", "FAIL"]
diffractionConditions = [
    soaksCollated["Diffraction Attempt %"] == 0,
    (soaksCollated["Diffraction Attempt %"] - soaksCollated["Diffraction Success %"]) >= diffractCutoffFailure,
]
soaksCollated["Diffraction Grade"] = np.select(diffractionConditions, diffractionChoices, default = "PASS")

numericResolution = pd.to_numeric(soaksCollated["Resolution (A)"], errors='coerce')  # avoids error at str to float comparison by replacing "NO DATA" with NaN in numericResolution
resolutionChoices = ["NO DATA", "FAIL", "PARTIAL"]
resolutionConditions = [
    soaksCollated["Resolution (A)"] == "NO DATA",
    numericResolution >= resolutionCutoffFailure,
    numericResolution >= resolutionCutoffPartial,
]
soaksCollated["Resolution Grade"] = np.select(resolutionConditions, resolutionChoices, default="PASS")

# Rcryst/DIMPLE Grade here - Fail cutoff is meant to find fail rates of 2-3 for 3 attempts, partial for 1 in 3. If there's lower diffraction success, will usually partial fail.
numericModelBuilding = soaksCollated[["Diffraction Success %", "DIMPLE Success %", "Rcryst Pass %", "Rfree Pass %"]].apply(pd.to_numeric, errors="coerce")
modelBuildingChoices = ["NO DATA", "FAIL", "FAIL", "FAIL", "PARTIAL", "PARTIAL", "PARTIAL"]
modelBuildingConditions = [
    soaksCollated["Diffraction Success %"] == 0,
    (numericModelBuilding["Diffraction Success %"] - numericModelBuilding["DIMPLE Success %"]) > modelBuildingCutoffFailure,
    (numericModelBuilding["Diffraction Success %"] - numericModelBuilding["Rcryst Pass %"]) > modelBuildingCutoffFailure,
    (numericModelBuilding["Diffraction Success %"] - numericModelBuilding["Rfree Pass %"]) > modelBuildingCutoffFailure,
    (numericModelBuilding["Diffraction Success %"] - numericModelBuilding["DIMPLE Success %"]) > modelBuildingCutoffPartial,
    (numericModelBuilding["Diffraction Success %"] - numericModelBuilding["Rcryst Pass %"]) > modelBuildingCutoffPartial,
    (numericModelBuilding["Diffraction Success %"] - numericModelBuilding["Rfree Pass %"]) > modelBuildingCutoffPartial,
]
soaksCollated["Model Building Grade"] = np.select(modelBuildingConditions, modelBuildingChoices, default="PASS")

gradeDf = soaksCollated[["Mount Grade", "Diffraction Grade",  "Resolution Grade",  "Model Building Grade"]]
def overallGrading(row):
    if any(row == "FAIL"):
        overallGrade = "FAIL"
    elif any(row == "PARTIAL"):
        overallGrade = "PARTIAL"
    elif all(row == "PASS"):
        overallGrade = "PASS"
    elif all(row == "NO DATA"):
        overallGrade = "NO DATA"
    else:
        overallGrade = "ERROR"
    return overallGrade

soaksCollated["Overall Grade"] = gradeDf.apply(overallGrading, axis=1)

# ---------------------- Improving sorting / Readability ----------------------

codeSplit = xceImport[["SoakCodes", "SolventFraction", "CompoundCode"]]
codeSplit = codeSplit.drop_duplicates(keep = "first")
soaksCollated = pd.merge(soaksCollated, codeSplit, how = "left", on = "SoakCodes", validate = "one_to_one")
soaksCollated = soaksCollated.rename(columns = {"CompoundCode": "Compound Code", "SolventFraction": "Solvent Fraction"})

def sortKey(col):
    if col.name == "Solvent Fraction":
        return natsort_keygen()(col)
    return col

soaksCollated = soaksCollated.sort_values(
    by = ["Compound Code", "Solvent Fraction"],
    key = sortKey,
    ascending=[True, False],
)

soaksCollated["In-Drop Molarity (mM)"] = (soaksCollated["Solvent Fraction"]/100)*compoundStockConcentration
commonTransferFrac = soaksCollated["Solvent Fraction"].mode()[0]
soaksCollated["Theoretical Source Molarity (mM)"] = ((soaksCollated["In-Drop Molarity (mM)"]) * (100 / commonTransferFrac))

# print(soaksCollated)
soaksCollated.to_csv("python_scripts/XtalControl/soakOutput.csv")

# ------------- table visualisation ----------------------
# most of this section is vibecoded

def formatDfText(inputvalue):
    if inputvalue == "NO DATA":
        return ("background-color: grey; color: white; text-align: center")
    elif inputvalue == "PASS":
        return ("background-color: green; color: white; text-align: center")
    elif inputvalue == "PARTIAL":
        return ("background-color: orange; color: white; text-align: center")
    elif inputvalue == "FAIL":
        return ("background-color: red; color: white; text-align: center")
    else: 
        return "" 

soaksCollated['_groupIndex'] = (soaksCollated['Compound Code'] != soaksCollated['Compound Code'].shift()).cumsum() % 2 # TRUE/FALSE row indicating when the Compound Code changes for formatting

def shadeGroups(row):
    if row['_groupIndex'] == 1:
        return ['background-color: #e8e8e8'] * len(row)
    return [''] * len(row)

(soaksCollated.style
    .format(precision=2, decimal=".")
    .apply(shadeGroups, axis=1)
    .map(formatDfText)
    .set_properties(**{
        'font-family': 'Arial, sans-serif',
        'font-size': '12px',
        'line-height': '24px',
        'white-space': 'nowrap',
        'border-bottom': '1px solid #ddd',
        'padding': '4px 8px',
    })
    .set_table_styles([
        {'selector': 'th', 'props': [
            ('font-family', 'Arial, sans-serif'),
            ('font-size', '12px'),
            ('padding', '4px 8px'),
            ('white-space', 'nowrap'),
            ('border-bottom', '2px solid #999'),
        ]},
        {'selector': 'th.col_heading', 'props': [
            ('background-color', 'white !important'),
            ('z-index', '3 !important'),
        ]},
    ], overwrite=False)
    .set_sticky(axis=1)
    .set_properties(subset=['Compound Code', 'Solvent Fraction', 'Theoretical Source Molarity (mM)'], **{'font-weight': 'bold'})
    .hide(subset=['_groupIndex'], axis='columns')
    .hide(axis="index")
    .to_html(buf="python_scripts/XtalControl/soakOutputHTML.html", bold_headers=True, doctype_html=True)
)

soaksCollated = soaksCollated.drop(columns='_groupIndex')