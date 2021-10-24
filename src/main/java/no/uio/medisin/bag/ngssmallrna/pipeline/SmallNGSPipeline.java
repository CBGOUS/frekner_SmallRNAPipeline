/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.StringWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import no.uio.medisin.bag.ngssmallrna.steps.NGSBase;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtieMapSingleReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepCollapseReads;
import no.uio.medisin.bag.ngssmallrna.steps.NGSRunStepJSONData;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStep;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStepSubclass;
import no.uio.medisin.bag.ngssmallrna.steps.StepParseSAMForMiRNAs;
import no.uio.medisin.bag.ngssmallrna.steps.StepDEwithEdgeR;
import no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeSAMforStartPositions;
import no.uio.medisin.bag.ngssmallrna.steps.StepBSMapCalcMethRatios;
import no.uio.medisin.bag.ngssmallrna.steps.StepBSMapReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtieMapPairedReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepDataCleanUp;
import no.uio.medisin.bag.ngssmallrna.steps.StepGrabFlankingSequence;
import no.uio.medisin.bag.ngssmallrna.steps.InputDataForStep;
import no.uio.medisin.bag.ngssmallrna.steps.NGSStepFactory;
import no.uio.medisin.bag.ngssmallrna.steps.StepAnalyzeSAMforCodingVsNonCoding;
import no.uio.medisin.bag.ngssmallrna.steps.StepBismarkMapReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepUnzipInputFiles;
import no.uio.medisin.bag.ngssmallrna.steps.StepSingleReadAdapterTrim;
import no.uio.medisin.bag.ngssmallrna.steps.StepBowtie2MapSingleReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepCalculateEntropyInBAM;
import no.uio.medisin.bag.ngssmallrna.steps.StepNGMMapSingleReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepConsolidatePreMiRNAData;
import no.uio.medisin.bag.ngssmallrna.steps.StepExit;
import no.uio.medisin.bag.ngssmallrna.steps.StepHiSatMapPairedReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepParseSAMForPreMiRNAReads;
import no.uio.medisin.bag.ngssmallrna.steps.StepFastqQC;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.yaml.snakeyaml.Yaml;
import org.codehaus.jackson.annotate.JsonAutoDetect;
import org.codehaus.jackson.annotate.JsonMethod;
import org.codehaus.jackson.map.ObjectMapper;
import org.yaml.snakeyaml.DumperOptions;

/**
 *
 * @author sr
 */
public class SmallNGSPipeline {

    
    static Logger                       logger           = LogManager.getLogger();
    static  final String                FILE_SEPARATOR   = System.getProperty("file.separator");
    static  final String                EXECUTE_BY_STEPS = "EXECUTE_BY_STEP";
    static  final String                EXECUTE_BY_DATA  = "EXECUTE_BY_DATA";
    static  final String                EXIT_PIPELINE    = "EXIT";
    
    private String                      configurationFile = "";
    private String                      pipelineFile = "";
    private String                      dataFile = "";
    static  Yaml                        yaml;
    
    private Boolean                     generateExampleConfigurationFile = false;
    
    private HashMap                     pipelineConfigurationDataHash;    
    private ReferenceDataLocations      refDataLocations;
    
    
    
    private ArrayList<String>           cleanupFiles    = new ArrayList<>();
    private PipelineData                pipelineData    = new PipelineData();
    private ArrayList<SampleDataEntry>  sampleData      = new ArrayList<>();

    private ArrayList<NGSBase>          ngsSteps        = new ArrayList<>();
    private String                      executeMode     = EXECUTE_BY_STEPS;
    
    
    /**
     * here we load the information from the three specified files
     * (run, pipeline and dataset) 
     * 
     * @throws IOException 
     */
    public void preparePipeline() throws IOException{      
        
        
        
        if(new File(this.getConfigurationFile()).exists()== false)
        {
            logger.error("configuration file <" + this.getConfigurationFile() + "> does not exist");
            throw new IOException("configuration file <" + this.getConfigurationFile() + "> does not exist");
        }

        if (new File(this.getPipelineFile()).exists() == false) {
            logger.error("pipeline file <" + this.getPipelineFile() + "> does not exist");
            throw new IOException("pipeline file <" + this.getPipelineFile() + "> does not exist");
        }

        if (new File(this.getDataFile()).exists() == false) {
            logger.error("data file <" + this.getDataFile() + "> does not exist");
            throw new IOException("data file <" + this.getDataFile() + "> does not exist");
        }

        this.readConfigurationFile();
        this.readPipelineFile();
        this.readDataFile();

    }

    /**
     *
     * This pre-processes the steps to ensure as far as possible that the
     * pipeline is runnable.
     *
     * @throws IOException, Exception
     */
    public void prepareSteps() throws IOException, Exception{
        
        if(this.executeMode.equals(EXECUTE_BY_DATA)){
            prepareStepsByData();
        }
        else{
            prepareStepsByStep();
        } 
    }

    /**
     * execute steps on the stack
     *
     * @throws IOException
     * @throws Exception
     */
    public void executePipeline() throws IOException, Exception {

        logger.info("executing pipeline by step list");
        for (NGSBase ngsStep : ngsSteps) {
            if(ngsStep == null) continue;
            ngsStep.verifyInputData();
            ngsStep.execute();
        }

    }

    /**
     * read file containing information about the data to be analyzed
     *
     * @throws IOException
     *
     */
    public void readDataFile() throws IOException {

        logger.info("\nread data file");
        String line = "";
        BufferedReader bwData = new BufferedReader(new FileReader(new File(this.getDataFile())));

        bwData.readLine(); // skip header line
        while ((line = bwData.readLine()) != null) {
            logger.info("-" + line);
            if (line.startsWith("#")) {
                continue;
            }

            String tokens[] = line.split("\t");
            String file1;
            String file2;
            if (tokens[0].contains(",")) {
                file1 = tokens[0].split(",")[0].trim();
                file2 = tokens[0].split(",")[1].trim();
            } else {
                file1 = tokens[0].split(",")[0].trim();
                file2 = null;
            }
            String source = tokens[1];
            String condition = tokens[2];
            String time = tokens[3];
            String note = "";
            if(line.split("\t").length == 5)
                note = tokens[4];

            getSampleData().add(new SampleDataEntry(file1, file2, source, condition, time, note));

        }
        bwData.close();
        logger.info("\nread " + getSampleData().size() + " entries");

    }

    /**
     * read configuration settings for the pipeline there should be settings for
     * each step. If the data is missing for a step the method should return a
     * warning.
     *
     * @throws IOException
     */
    public void readConfigurationFile() throws IOException {

        logger.info("read pipeline configuration file <" + getConfigurationFile() + ">");

        yaml = new Yaml();
        pipelineConfigurationDataHash = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));

        HashMap refDataOptions = (HashMap) pipelineConfigurationDataHash.get(ReferenceDataLocations.ID_CONFIG_ID);
        if (refDataOptions == null) {
            logger.info("yaml configuration file is missing reference data locations");
            throw new IOException("yaml configuration file is missing reference data locations");
        }
        refDataLocations = new ReferenceDataLocations(refDataOptions);
        if (refDataLocations == null) {
            logger.info("failed to load the reference data locations correctly");
            throw new IOException("failed to load the reference data locations correctly");
        }

        logger.info("done\n");

    }


    /**
     * read configuration settings for the pipeline from XML format
     * there should be settings for each step
     * If the data is missing for a step the method should return a warning
     *
     * @throws IOException
     */
    public void readXMLConfigurationFile() throws IOException {

        logger.info("read pipeline XML configuration file <" + getConfigurationFile() + ">");

        yaml = new Yaml();
        pipelineConfigurationDataHash = (HashMap) yaml.load(new FileInputStream(new File(getConfigurationFile())));

        HashMap refDataOptions = (HashMap) pipelineConfigurationDataHash.get(ReferenceDataLocations.ID_CONFIG_ID);
        refDataLocations = new ReferenceDataLocations(refDataOptions);
        if (refDataLocations == null) {
            logger.info("failed to load the reference data locations correctly");
            throw new IOException("failed to load the reference data locations correctly");
        }

        logger.info("done\n");

    }


    /**
     * generate a sample configuration file. this is primarily to inform users
     * what can be specified for each step
     *
     * @throws IOException
     */
    public void generateSampleConfigurationFile() throws IOException {

        logger.info("write example pipeline configuration file");
        DumperOptions options = new DumperOptions();
        options.setDefaultFlowStyle(DumperOptions.FlowStyle.BLOCK);

        yaml = new Yaml(options);

        InputDataForStep emptySID = new InputDataForStep("", "", new ReferenceDataLocations(), "", "", "", new ArrayList<SampleDataEntry>());
        Map<String, Object> pipelineExampleConfiguration = new HashMap();

        ReferenceDataLocations rdl = new ReferenceDataLocations();
        pipelineExampleConfiguration.put(ReferenceDataLocations.ID_CONFIG_ID, rdl.generateExampleConfigurationData());

        StepUnzipInputFiles stepUnzip = new StepUnzipInputFiles(emptySID);
        pipelineExampleConfiguration.put(StepUnzipInputFiles.STEP_ID_STRING, stepUnzip.generateExampleConfigurationData());

        StepSingleReadAdapterTrim stepSingleAdapterTrim = new StepSingleReadAdapterTrim(emptySID);
        pipelineExampleConfiguration.put(StepSingleReadAdapterTrim.STEP_ID_STRING, stepSingleAdapterTrim.generateExampleConfigurationData());

        StepCollapseReads stepcollapseReads = new StepCollapseReads(emptySID);
        pipelineExampleConfiguration.put(StepCollapseReads.STEP_ID_STRING, stepcollapseReads.generateExampleConfigurationData());

        StepBowtieMapSingleReads stepBowtieSingleMap = new StepBowtieMapSingleReads(emptySID);
        pipelineExampleConfiguration.put(StepBowtieMapSingleReads.STEP_ID_STRING, stepBowtieSingleMap.generateExampleConfigurationData());
        
        StepBowtie2MapSingleReads stepBowtie2SingleMap = new StepBowtie2MapSingleReads(emptySID);
        pipelineExampleConfiguration.put(StepBowtie2MapSingleReads.STEP_ID_STRING, stepBowtie2SingleMap.generateExampleConfigurationData());
        //pipelineExampleConfiguration.put(StepBowtie2MapSingleReads.STEP_ID_STRING, StepBowtie2MapSingleReads.generateExampleConfigurationData());

        StepHiSatMapPairedReads stepHiSatPairedMap = new StepHiSatMapPairedReads(emptySID);
        pipelineExampleConfiguration.put(StepHiSatMapPairedReads.STEP_ID_STRING, stepHiSatPairedMap.generateExampleConfigurationData());

        StepNGMMapSingleReads stepNGMSingleMap = new StepNGMMapSingleReads(emptySID);
        pipelineExampleConfiguration.put(StepNGMMapSingleReads.STEP_ID_STRING, stepNGMSingleMap.generateExampleConfigurationData());

        StepCalculateEntropyInBAM stepCalculateEntropyInBAM = new StepCalculateEntropyInBAM(emptySID);
        pipelineExampleConfiguration.put(StepCalculateEntropyInBAM.STEP_ID_STRING, stepCalculateEntropyInBAM.generateExampleConfigurationData());

        StepBSMapReads stepBSMapReads = new StepBSMapReads(emptySID);
        pipelineExampleConfiguration.put(StepBSMapReads.STEP_ID_STRING, stepBSMapReads.generateExampleConfigurationData());

        StepBismarkMapReads stepBismarkMapReads = new StepBismarkMapReads(emptySID);
        pipelineExampleConfiguration.put(StepBismarkMapReads.STEP_ID_STRING, stepBismarkMapReads.generateExampleConfigurationData());

        StepBSMapCalcMethRatios stepCalcMethRatios = new StepBSMapCalcMethRatios(emptySID);
        pipelineExampleConfiguration.put(StepBSMapCalcMethRatios.STEP_ID_STRING, stepCalcMethRatios.generateExampleConfigurationData());

        StepAnalyzeSAMforStartPositions stepAnalyzeStartPos = new StepAnalyzeSAMforStartPositions(emptySID);
        pipelineExampleConfiguration.put(StepAnalyzeSAMforStartPositions.STEP_ID_STRING, stepAnalyzeStartPos.generateExampleConfigurationData());

        StepAnalyzeSAMforCodingVsNonCoding stepAnalyzeStartPosCodingVsNonCoding = new StepAnalyzeSAMforCodingVsNonCoding(emptySID);
        pipelineExampleConfiguration.put(StepAnalyzeSAMforCodingVsNonCoding.STEP_ID_STRING, stepAnalyzeStartPosCodingVsNonCoding.generateExampleConfigurationData());

        StepParseSAMForMiRNAs stepParseSAMforMirs = new StepParseSAMForMiRNAs(emptySID);
        pipelineExampleConfiguration.put(StepParseSAMForMiRNAs.STEP_ID_STRING, stepParseSAMforMirs.generateExampleConfigurationData());

        StepParseSAMForPreMiRNAReads stepParseSAMforPreMirReads = new StepParseSAMForPreMiRNAReads(emptySID);
        pipelineExampleConfiguration.put(StepParseSAMForPreMiRNAReads.STEP_ID_STRING, stepParseSAMforPreMirReads.generateExampleConfigurationData());

        StepConsolidatePreMiRNAData stepConsolidatePreMiRNAData = new StepConsolidatePreMiRNAData(emptySID);
        pipelineExampleConfiguration.put(StepConsolidatePreMiRNAData.STEP_ID_STRING, stepConsolidatePreMiRNAData.generateExampleConfigurationData());

        StepDEwithEdgeR stepDEwithEdgeR = new StepDEwithEdgeR(emptySID);
        pipelineExampleConfiguration.put(StepDEwithEdgeR.STEP_ID_STRING, stepDEwithEdgeR.generateExampleConfigurationData());

        StepGrabFlankingSequence stepGrabFlankingSequence = new StepGrabFlankingSequence(emptySID);
        pipelineExampleConfiguration.put(StepGrabFlankingSequence.STEP_ID_STRING, stepGrabFlankingSequence.generateExampleConfigurationData());
        
        StepBowtieMapPairedReads stepBowtieMapPairedReads = new StepBowtieMapPairedReads(emptySID);
        pipelineExampleConfiguration.put(StepBowtieMapPairedReads.STEP_ID_STRING, stepBowtieMapPairedReads.generateExampleConfigurationData());

        StepDataCleanUp stepCleanUp = new StepDataCleanUp(emptySID);
        pipelineExampleConfiguration.put(StepDataCleanUp.STEP_ID_STRING, stepCleanUp.generateExampleConfigurationData());
        
        StepFastqQC stepFastQC = new StepFastqQC(emptySID);
        pipelineExampleConfiguration.put(StepFastqQC.STEP_ID_STRING, stepFastQC.generateExampleConfigurationData());

        //String dumpFilename = System.getProperty("user.dir") + FILE_SEPARATOR + "pipelineConfiguration.sample.yaml";
        String dumpFilename = this.getConfigurationFile();
        StringWriter writer = new StringWriter();
        yaml.dump(pipelineExampleConfiguration, writer);
        logger.info(writer);
        try (FileWriter sampleFileWriter = new FileWriter(new File(dumpFilename))) {
            yaml.dump(pipelineExampleConfiguration, sampleFileWriter);
        }
    }

    /**
     * read pipeline file
     *
     * @throws IOException
     */
    public void readPipelineFile() throws IOException {

        ObjectMapper mapper = new ObjectMapper();
        mapper.setVisibility(JsonMethod.FIELD, JsonAutoDetect.Visibility.ANY);

        logger.info("read pipeline file <" + this.getPipelineFile() + ">");
        try {
            pipelineData = mapper.readValue(new File(this.getPipelineFile()), PipelineData.class);
            if(!pipelineData.isDataComplete()){
                logger.error("pipeline data has missing fields:");
                logger.error(pipelineData.toString());
                throw new RuntimeException("pipeline data has missing fields:\n" + pipelineData.toString());
            }            
        } catch (IOException ex) {
            logger.error("error reading <" + this.getPipelineFile() + ">");
            logger.error(ex.toString());
            throw new IOException("error reading <" + this.getPipelineFile() + ">");
        }
        logger.info("done");
        logger.info(getPipelineData().toString());

    }

    /**
     * @return the configurationFile
     */
    public String getConfigurationFile() {
        return configurationFile;
    }

    /**
     * @param configurationFile the configurationFile to set
     */
    public void setConfigurationFile(String configurationFile) {
        this.configurationFile = configurationFile;
    }

    /**
     * @return the pipelineFile
     */
    public String getPipelineFile() {
        return pipelineFile;
    }

    /**
     * @param pipelineFile the pipelineFile to set
     */
    public void setPipelineFile(String pipelineFile) {
        this.pipelineFile = pipelineFile;
    }

    /**
     * @return the dataFile
     */
    public String getDataFile() {
        return dataFile;
    }

    /**
     * @param dataFile the dataFile to set
     */
    public void setDataFile(String dataFile) {
        this.dataFile = dataFile;
    }

    /**
     * @return the pipelineData
     */
    public PipelineData getPipelineData() {
        return pipelineData;
    }

    /**
     * @return the SampleData
     */
    public ArrayList<SampleDataEntry> getSampleData() {
        return sampleData;
    }

    /**
     * @return the cleanupFiles
     */
    public ArrayList<String> getCleanupFiles() {
        return cleanupFiles;
    }

    /**
     * @param cleanupFiles the cleanupFiles to set
     */
    public void setCleanupFiles(ArrayList<String> cleanupFiles) {
        this.cleanupFiles = cleanupFiles;
    }

    /**
     * @return the dataLocations
     */
    public ReferenceDataLocations getDataLocations() {
        return refDataLocations;
    }

    /**
     * @param dataLocations the dataLocations to set
     */
    public void setDataLocations(ReferenceDataLocations dataLocations) {
        this.refDataLocations = dataLocations;
    }

    /**
     * @return the generateExampleConfigurationFile
     */
    public Boolean getGenerateExampleConfigurationFile() {
        return generateExampleConfigurationFile;
    }

    /**
     * @param generateExampleConfigurationFile the
     * generateExampleConfigurationFile to set
     */
    public void setGenerateExampleConfigurationFile(Boolean generateExampleConfigurationFile) {
        this.generateExampleConfigurationFile = generateExampleConfigurationFile;
    }

//<<<<<<< HEAD
    /**
     * set execute model defined by user,  execute by data or by steps
     * @param executeModelString set
     */
    public void setExecuteModel(String executeModelString) {
        this.executeMode = executeModelString;
    }

    
    /**
     * The default approach is to execute by Step, i.e. all data is pushed 
     * through each step at one time. However, in some circumstances, the user
     * may prefer to execute by data - i.e., all the steps are executed for 
     * a single data file, then the loop moves to the next data file. 
     * This can help the pipeline run within a smaller disk footprint.
     * (e.g., when analyzing paired-end mRNA data, large quantities
     * of disk space will be consumed and only released when the cleanup step
     * is executed). However, this makes no logical sense for steps such as
     * Different Expression analysis or Merging Count Files which cannot be
     * executed until all the required input data has been generated.
     * 
     * 
     * @throws Exception 
     */
    private void prepareStepsByData() throws Exception {
        logger.info("preparing steps for datawise execution");
                
        List <NGSRunStepJSONData> unDatableStepsData = new ArrayList<>();
        logger.info("adding DATABLE steps ");
        
        int entryCount = 1;
        for(SampleDataEntry thisSampleDataEntry: sampleData){
            logger.info(thisSampleDataEntry.getFastqFile1());
            for (NGSRunStepJSONData stepData: this.getPipelineData().getStepsJSONData()){
                logger.info(" - Found step " + stepData.getStep());
                if(stepData.getStep().toUpperCase().equals(StepExit.STEP_ID_STRING))
                    break;
                NGSStepSubclass thisStepSubclass 
                        =  ((NGSStep)Class.forName("no.uio.medisin.bag.ngssmallrna.steps.Step" 
                                + stepData.getStep()).newInstance()).classSubtype;
                
                if(thisStepSubclass == NGSStepSubclass.DATABLE){
                    logger.info(" - step is DATABLE");
                    ArrayList<SampleDataEntry> thisSampleData = new ArrayList<>();
                    thisSampleData.add(thisSampleDataEntry);
                    InputDataForStep stepInputData = new InputDataForStep(this.getPipelineData().getProjectID(), 
                            this.getPipelineData().getProjectRoot(), 
                            refDataLocations, 
                            stepData.getInputFolder(), 
                            stepData.getOutputFolder(), 
                            stepData.getParameterString(),
                            this.getSampleData());   

                    logger.info("loading step " + stepData.getStep());

                    NGSStep newStep = NGSStepFactory.getNewStep(stepData.getStep());

                    newStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(stepData.getStep()));
                    newStep.setStepInputData(stepInputData);
                    ngsSteps.add(newStep);
                }else{
                    logger.info(" - step is UNDATABLE, delay");
                    if(entryCount==1){ // only need to add an UNDATABLE step once
                        unDatableStepsData.add(stepData);                        
                    }
                }
                
            }
            entryCount++;
            
        }
        
        // add the steps that are undatable
        logger.info("adding UNDATABLE steps ");
        for (NGSRunStepJSONData stepData: unDatableStepsData){
            logger.info(" - Found step " + stepData.getStep());
            if(stepData.getStep().equals(StepExit.STEP_ID_STRING))
                return;
            InputDataForStep stepInputData = new InputDataForStep(this.getPipelineData().getProjectID(), 
                    this.getPipelineData().getProjectRoot(), 
                    refDataLocations, 
                    stepData.getInputFolder(), 
                    stepData.getOutputFolder(), 
                    stepData.getParameterString(),
                    this.getSampleData());   
            logger.info("loading step " + stepData.getStep());
            
            NGSStep newStep = NGSStepFactory.getNewStep(stepData.getStep());
            
            newStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(stepData.getStep()));
            newStep.setStepInputData(stepInputData);
            ngsSteps.add(newStep);
        }
    
        logger.info("completed");
    }

    
    
    
    /**
     * This executes the pipeline one step at a time and loops through all
     * input data for that step. The drawback is that large amounts of disk space
     * can be used for steps such as mapping and the space is only released
     * at the end of the pipeline when the Cleanup step is executed
     * 
     * @throws Exception 
     */
    private void prepareStepsByStep() throws Exception {
        logger.info("preparing steps for stepwise execution");
        for (NGSRunStepJSONData stepJSONData: this.getPipelineData().getStepsJSONData()){
            logger.info(" - Found step " + stepJSONData.getStep());
            if(stepJSONData.getStep().equals(StepExit.STEP_ID_STRING))
                return;
            InputDataForStep stepInputData = new InputDataForStep(this.getPipelineData().getProjectID(), 
                    this.getPipelineData().getProjectRoot(), 
                    refDataLocations, 
                    stepJSONData.getInputFolder(), 
                    stepJSONData.getOutputFolder(), 
                    stepJSONData.getParameterString(),
                    this.getSampleData());   
            logger.info("loading step " + stepJSONData.getStep());
            
            NGSStep newStep = NGSStepFactory.getNewStep(stepJSONData.getStep());
            
            newStep.parseConfigurationData((HashMap) pipelineConfigurationDataHash.get(stepJSONData.getStep()));
            newStep.setStepInputData(stepInputData);
            newStep.parseStepParameters();
            ngsSteps.add(newStep);
        }
    }

}    
    
