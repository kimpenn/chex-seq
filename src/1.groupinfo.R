source("src/lib/CHEX.R")
source("src/lib/Repo.R")
source("src/lib/NGS.R")

SampleInfo <- read.csv("data/SampleInfo.csv", as.is = TRUE, check.names = FALSE)
sampleIDs <- SampleInfo$SampleID
cat(sampleIDs, file = "data/sampleIDs.txt", sep = "\n")

bioGroups <- c(
    "K562", "K562TPAnone", "K562TPA15min", "K562TPA1hr", "K562TPA2hr", "K562TPA24hr", 
    "HumanAstroCulture", "HumanNeuronCulture", "HumanInterneuronCulture",
    "MouseAstroCulture", "MouseNeuronCulture", "MouseNeuronSlice", "MouseInterneuronSlice", 
    "MouseLungEpithelialCulture", "MouseLungEndothelialCulture",
    "MouseKidneyEpithelialCulture",
    "MouseCardiomyoCulture", "MouseCardiomyoSlice",
    "RatCardiomyoCulture", 
    "NoCell", 
    "K562MungBean", 
    "HBR"
)
bioGroup2species = c(
    K562 = "human", K562TPAnone = "human", K562TPA15min = "human", K562TPA1hr = "human", K562TPA2hr = "human", K562TPA24hr = "human",
    HumanAstroCulture = "human", HumanNeuronCulture = "human", HumanInterneuronCulture = "human",
    MouseAstroCulture = "mouse", MouseNeuronCulture = "mouse", MouseNeuronSlice = "mouse", MouseInterneuronSlice = "mouse",
    MouseLungEpithelialCulture = "mouse", MouseLungEndothelialCulture = "mouse", 
    MouseKidneyEpithelialCulture = "mouse",
    MouseCardiomyoCulture = "mouse", MouseCardiomyoSlice = "mouse", 
    RatCardiomyoCulture = "rat", 
    NoCell = "none",
    K562MungBean = "human",
    HBR = "human"
)

bioGroup2cellType = c(
    K562 = "K562", K562TPAnone = "K562", K562TPA15min = "K562", K562TPA1hr = "K562", K562TPA2hr = "K562", K562TPA24hr = "K562",
    HumanAstroCulture = "Astro", HumanNeuronCulture = "Neuron", HumanInterneuronCulture = "Interneuron",
    MouseAstroCulture = "Astro", MouseNeuronCulture = "Neuron", MouseNeuronSlice = "Neuron", MouseInterneuronSlice = "Interneuron",
    MouseLungEpithelialCulture = "LungEpithelial", MouseLungEndothelialCulture = "LungEndothelial",
    MouseKidneyEpithelialCulture = "KidneyEpithelial",
    MouseCardiomyoCulture = "Cardiomyo", MouseCardiomyoSlice = "Cardiomyo",
    RatCardiomyoCulture = "Cardiomyo",
    NoCell = "none",
    K562MungBean = "K562",
    HBR = "Brain"
)

bioGroup2harvest = c(
    K562 = "Culture", K562TPAnone = "Culture", K562TPA15min = "Culture", K562TPA1hr = "Culture", K562TPA2hr = "Culture", K562TPA24hr = "Culture",
    HumanAstroCulture = "Culture", HumanNeuronCulture = "Culture", HumanInterneuronCulture = "Culture",
    MouseAstroCulture = "Culture", MouseNeuronCulture = "Culture", MouseNeuronSlice = "Slice", MouseInterneuronSlice = "Slice",
    MouseLungEpithelialCulture = "Culture", MouseLungEndothelialCulture = "Culture", 
    MouseKidneyEpithelialCulture = "Culture",
    MouseCardiomyoCulture = "Culture", MouseCardiomyoSlice = "Slice",
    RatCardiomyoCulture = "Culture", 
    NoCell = "none",
    K562MungBean = "Culture",
    HBR = "none"
)

bioGroup2drugGroup = c(
    K562 = "none", K562TPAnone = "TPA", K562TPA15min = "TPA", K562TPA1hr = "TPA", K562TPA2hr = "TPA", K562TPA24hr = "TPA",
    HumanAstroCulture = "none", HumanNeuronCulture = "none", HumanInterneuronCulture = "none",
    MouseAstroCulture = "none", MouseNeuronCulture = "none", MouseNeuronSlice = "none", MouseInterneuronSlice = "none",
    MouseLungEpithelialCulture = "none", MouseLungEndothelialCulture = "none",
    MouseKidneyEpithelialCulture = "none",
    MouseCardiomyoCulture = "none", MouseCardiomyoSlice = "none", 
    RatCardiomyoCulture = "none", 
    NoCell = "none",
    K562MungBean = "none",
    HBR = "none"
)

bioGroup2hasTPA = c(
    K562 = "N", K562TPAnone = "N", K562TPA15min = "Y", K562TPA1hr = "Y", K562TPA2hr = "Y", K562TPA24hr = "Y",
    HumanAstroCulture = "N", HumanNeuronCulture = "N", HumanInterneuronCulture = "N",
    MouseAstroCulture = "N", MouseNeuronCulture = "N", MouseNeuronSlice = "N", MouseInterneuronSlice = "N",
    MouseLungEpithelialCulture = "N", MouseLungEndothelialCulture = "N",
    MouseKidneyEpithelialCulture = "N",
    MouseCardiomyoCulture = "N", MouseCardiomyoSlice = "N", 
    RatCardiomyoCulture = "N", 
    NoCell = "N",
    K562MungBean = "N",
    HBR = "N"
)

bioGroup2TPADose = c(
    K562 = "none", K562TPAnone = "none", K562TPA15min = "16uM", K562TPA1hr = "16uM", K562TPA2hr = "16uM", K562TPA24hr = "16uM",
    HumanAstroCulture = "none", HumanNeuronCulture = "none", HumanInterneuronCulture = "none", 
    MouseAstroCulture = "none", MouseNeuronCulture = "none", MouseNeuronSlice = "none", MouseInterneuronSlice = "none",
    MouseLungEpithelialCulture = "none", MouseLungEndothelialCulture = "none",
    MouseKidneyEpithelialCulture = "none",
    MouseCardiomyoCulture = "none", MouseCardiomyoSlice = "none", 
    RatCardiomyoCulture = "none", 
    NoCell = "none",
    K562MungBean = "none",
    HBR = "none"
)

bioGroup2TPATime = c(
    K562 = "none", K562TPAnone = "none", K562TPA15min = "15min", K562TPA1hr = "1hr", K562TPA2hr = "2hr", K562TPA24hr = "24hr",
    HumanAstroCulture = "none", HumanNeuronCulture = "none", HumanInterneuronCulture = "none", 
    MouseAstroCulture = "none", MouseNeuronCulture = "none", MouseNeuronSlice = "none", MouseInterneuronSlice = "none",
    MouseLungEpithelialCulture = "none", MouseLungEndothelialCulture = "none",
    MouseKidneyEpithelialCulture = "none",
    MouseCardiomyoCulture = "none", MouseCardiomyoSlice = "none", 
    RatCardiomyoCulture = "none", 
    NoCell = "none",
    K562MungBean = "none",
    HBR = "none"
)

bioGroup2hasMungBean = c(
    K562 = "N", K562TPAnone = "N", K562TPA15min = "N", K562TPA1hr = "N", K562TPA2hr = "N", K562TPA24hr = "N",
    HumanAstroCulture = "N", HumanNeuronCulture = "N", HumanInterneuronCulture = "N", 
    MouseAstroCulture = "N", MouseNeuronCulture = "N", MouseNeuronSlice = "N", MouseInterneuronSlice = "N",
    MouseLungEpithelialCulture = "N", MouseLungEndothelialCulture = "N",
    MouseKidneyEpithelialCulture = "N",
    MouseCardiomyoCulture = "N", MouseCardiomyoSlice = "N", 
    RatCardiomyoCulture = "N", 
    NoCell = "N",
    K562MungBean = "Y",
    HBR = "N"
)

CHEX$createBioGroupConfig(
    bioGroups = bioGroups, 
    bioGroup2species = bioGroup2species, 
    bioGroup2cellType = bioGroup2cellType, 
    bioGroup2harvest = bioGroup2harvest, 
    bioGroup2drugGroup = bioGroup2drugGroup, 
    bioGroup2hasTPA = bioGroup2hasTPA, 
    bioGroup2TPADose = bioGroup2TPADose, 
    bioGroup2TPATime = bioGroup2TPATime, 
    bioGroup2hasMungBean = bioGroup2hasMungBean, 
    outFile = "data/BioGroupConfig.csv"
)

## E.639, 640, 641 are done on "mid-throughput" flowcells, which are resequenced in E.642, 648, 649.
## We don't need to include the mid-throughput data as they are technically replicates and hence add no biological information.
NoUseExptIDs <- c(639, 640, 641, 672)
## Aggregate single cells, bulk samples of the same biological group and the same experiment condition into 'virtual' samples
CHEX$createSampleInfoFull(
    sampleInfoFile = "data/SampleInfo.csv", 
    sampleInfoFullFile = "data/SampleInfoFull.csv",
    bioGroups = bioGroups,
    bioGroupConfigFile = "data/BioGroupConfig.csv",  
    sampleIDsOutFile = NULL,
    exptIDsToExclude = NoUseExptIDs
)
