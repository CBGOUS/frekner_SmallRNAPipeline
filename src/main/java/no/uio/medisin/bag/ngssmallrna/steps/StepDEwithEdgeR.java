/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.security.SecureRandom;
import java.util.Date;
import java.util.HashMap;
import java.util.Random;

import no.uio.medisin.bag.core.mirna.MiRNA;
import no.uio.medisin.bag.core.mirna.MiRNASet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import static no.uio.medisin.bag.ngssmallrna.steps.StepParseSAMForMiRNAs.logger;
import org.apache.commons.lang3.StringUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   Perform differential expression analysis using EdgeR

 * 
 * @author sr
 */

public class StepDEwithEdgeR extends NGSStep{
    
    static Logger                       logger                      = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.NOTDATABLE;

    public  static final String     STEP_ID_STRING                  = "DEwithEdgeR";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_PVALUE                       = "pValue";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";

    private static final String     GROUPS_FILE_EXTENSION           = ".groups.tsv";
    private static final String     MIR_COUNTS_EXTENSION            = ".trim.clp.gen.mircounts.tsv";
    private static final String     DE_RESULTS_EXTENSION            = ".DE.edgeR.sort.csv";
    private static final String     DE_SUMMARY_COUNTS_EXTENSION     = ".DE.edgeR.cbs.csv";
    private static final String     DE_SUMMARY_EXTENSION            = ".DE.edgeR.summary.txt";
    private static final String     DE_MDS_PLOT_EXTENSION           = ".DE.edgeR.MDS.";
    private static final String     DE_BCV_PLOT_EXTENSION           = ".DE.edgeR.BCV.";
    private static final String     DE_SMEAR_PLOT_EXTENSION         = ".DE.edgeR.Smear.";
    
    private static final String     RAW_INPUT_EXTENSION             = ".fastq.gz";
    
    private String                  miRBaseRelease                  = "22.1";
    private String                  ReferenceGenome                 = "";
    private double                  pValue                          = 0.05;
        
    private              String     rScriptFilename                 = "";
        
    private String                  mergedCountsFile;
    private String                  groupsFile;
    private String                  deResultsFile;
    private String                  deCountsBySampleFile;
    private String                  deSummaryFile;
    private String                  dePlotBCVfile;
    private String                  dePlotMDSfile;
    private String                  dePlotSmearfile;
    
    private static final int        PLOT_WIDTH                      = 8;
    private static final int        PLOT_HEIGHT                     = 8;
    private static final String     PLOT_UNITS                      = "in";
    private static final int        PLOT_RESOLUTION                 = 300;
    
    
    private MiRNASet           miRBaseMiRNAList                = new MiRNASet();

    
    
    public StepDEwithEdgeR(){
        classSubtype = NGSStepSubclass.NOTDATABLE;
    }
    
    public StepDEwithEdgeR(InputDataForStep sid){
        classSubtype = NGSStepSubclass.NOTDATABLE;
        stepInputData = sid;
    }



    @Override
    public String shortStepDescription(){
      return "Perform differential expression analysis using EdgeR";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "Perform differential expression analysis using EdgeR.\n"
              + "requires count file as input.\n";
    }


    @Override
    public void parseStepParameters() throws Exception{
      
    }
    


    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
    private static final String     ID_REF_GENOME                   = "host:";
    private static final String     ID_PVALUE                       = "pValue:";
    private static final String     ID_MIRBASE_VERSION              = "mirbase_release:";
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{
        
        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_PVALUE)==null) {
            logger.info("<" + ID_PVALUE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_PVALUE + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }

        
        try{
            this.setpValue((Double) configData.get(ID_PVALUE));
        }
        catch(NumberFormatException exNm){
            logger.info(ID_PVALUE + " <" + configData.get(ID_PVALUE) + "> is not an integer");
            throw new NumberFormatException(ID_PVALUE + " <" + configData.get(ID_PVALUE) + "> is not an integer");
        }        
        if (this.getpValue() <= 0 || this.getpValue() > 1.0){
            logger.info(ID_PVALUE + " <" + configData.get(ID_PVALUE) + "> must be an float between 0.0 and 1.0");
            throw new IllegalArgumentException(ID_PVALUE + " <" + configData.get(ID_PVALUE) + "> must be an float between 0.0 and 1.0");
        }
              
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        
        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }
        

        logger.info("passed");
    }
    
    
    
    @Override
    public void execute() throws IOException{
        
        logger.info(STEP_ID_STRING + ": execute step");        
        
        /*
            1. read in all sample count files and merge
            2. output merged count file
            3. generate R script to perform DE using EdgeR
            4. process EdgeR output file 
        */
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        
//        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
//                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
//        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
//        miRBaseMiRNAList.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);

//        logger.info("Merging Count Files");
//        String headerLine = "name";
//        String[] countStrings = new String[miRBaseMiRNAList.getNumberOfEntries()];
//        Arrays.fill(countStrings, "");
//        
//        Iterator itSD = this.getStepInputData().getSampleData().iterator();
//        while (itSD.hasNext()){
//            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
//            headerLine = headerLine.concat("\t" + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ""));
//            String  miRCountsFile  = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MIR_COUNTS_EXTENSION));
//            try{
//                int m=0;
//                BufferedReader brmiRCounts  = new BufferedReader(new FileReader(new File(miRCountsFile)));
//                    String countLine = "";
//                    while((countLine = brmiRCounts.readLine()) != null){
//                        countStrings[m] = countStrings[m].concat("\t" + countLine.split("\t")[1].trim() );
//                        m++;
//                    }
//                brmiRCounts.close();
//            }
//            catch(IOException ex){
//                logger.error("error reading count files for merging <" + miRCountsFile + "> \n" + ex.toString());
//                throw new IOException(STEP_ID_STRING + " :error reading counts file  <" + miRCountsFile + ">");
//            }
//        }
        
//        logger.info("Writing merged count files");
        //need add verify code in 
        mergedCountsFile        = inFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".merged.mirna_counts.tsv"; 
//        
//        try{
//            BufferedWriter bwMC = new BufferedWriter(new FileWriter(new File(mergedCountsFile)));
//            
//            bwMC.write(headerLine + "\n");
//            int m=0;
//            for(MiRNA miR: miRBaseMiRNAList.getMiRBaseMiRNAList()){
//                bwMC.write(miR.getMimatID() + "|" + miR.getName() + countStrings[m] + "\n");
//                m++;
//            }
//            
//            bwMC.close();
//        }
//        catch(IOException exIO){
//            logger.info("error writing merged counts File <" + mergedCountsFile + ">\n" + exIO);        
//            throw new IOException(STEP_ID_STRING + " :error writing merged counts file  <" + mergedCountsFile + ">");
//        }
        generateGroupsFile();
        buildRScript();
        executeRScript();
        
        
        logger.info(STEP_ID_STRING + ": completed");
        
    }
    
    
    
    /**
     * 
     * generates the a file that contains all the grouping information for all
     * samples in the experiment.
     * 
     * Grouping is according to the Condition column in the data file

          this has the format:
            File               Source	Condition   Time	Note
            SRR1642941.fastq	P1	U           NA	(U|44|CR|M|IF)
            SRR1642942.fastq	P1	T           NA	(T|44|CR|M|IF)
            SRR1642943.fastq	P2	U           NA	(U|52|NC|M|IF)
            SRR1642944.fastq	P2	T           NA	(T|51|NC|M|IF)
        
            i.e., Condition is equivalent to Group
        
            Need the file in the format
            Group	U   T   U   T
            sample names    SRR1642941  SRR1642942  SRR1642943  SRR1642944
            
            We could write the R to parse the sample file, but it makes the code
            harder to understand

     * 
     */
    private void generateGroupsFile() throws IOException{
        

        String groupString = "Group";
        String sampleString = "sample names";

        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            groupString = groupString.concat("\t" + sampleData.getCondition());
            sampleString = sampleString.concat("\t" + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ""));
        }        
        
        groupsFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + GROUPS_FILE_EXTENSION;
        logger.info("writing groups file "  + groupsFile);
        try{
            BufferedWriter bwGF = new BufferedWriter(new FileWriter(new File(groupsFile)));    
                bwGF.write(groupString + "\n");
                bwGF.write(sampleString + "\n");
            bwGF.close();
        }
        catch(IOException exGF){
            logger.error("error writing groups file "  + groupsFile);
            logger.error(exGF);
            throw new IOException(STEP_ID_STRING + " :error writing groups file "  + groupsFile);
        }

                
    }
    
    
    
    
    
    /**
     * construct the R script to perform DE analysis for this dataset
     * Need to create a new R script for each analysis and give it a
     * unique name so we can go back and check what parameters were
     * used in the analysis.
     * 
     * need to generate the input files based on the specified groups
     * 
     * @throws IOException
     */
    private void buildRScript() throws IOException{


        BigInteger big = new BigInteger(130, new Random());
        String randomName = new BigInteger(130, new SecureRandom()).toString(32);
        rScriptFilename = outFolder + FILESEPARATOR + randomName + ".R";
        rScriptFilename = rScriptFilename.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
        
        deResultsFile           = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + DE_RESULTS_EXTENSION;    
        deCountsBySampleFile    = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + DE_SUMMARY_COUNTS_EXTENSION;
        deSummaryFile           = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + DE_SUMMARY_EXTENSION;
        
        int minCounts = 10;
        ArrayList<String> cmdSet = new ArrayList<>();

        cmdSet.add("#");
        cmdSet.add("#   Rscript generated for project " + this.getStepInputData().getProjectID());
        cmdSet.add("#   created: " + new Date());
        cmdSet.add("#");
        cmdSet.add("#");
        cmdSet.add("library(edgeR)");
        cmdSet.add("RawCounts<-read.delim(\"" + mergedCountsFile + "\", row.names=\"name\")");
        cmdSet.add("g<-read.csv(\"" + groupsFile + "\", header=FALSE, sep=\"\\t\", row.names=1)");
        cmdSet.add("sampleGroups <- unname(unlist(g[\"Group\",]))");
        cmdSet.add("CountsDGE <- DGEList(counts=RawCounts, group=sampleGroups)");
        cmdSet.add("");
        cmdSet.add("CountsAboveNoise <- rowSums(cpm(CountsDGE)>" + minCounts + ") >= 2");
        cmdSet.add("CountsLessNoiseDGE <- CountsDGE[CountsAboveNoise,]");
        cmdSet.add("CountsLessNoiseDGE$samples$lib.size <- colSums(CountsLessNoiseDGE$counts)");
        cmdSet.add("");
        cmdSet.add("NormFactors <- calcNormFactors(CountsLessNoiseDGE)");
        cmdSet.add("");
        cmdSet.add("CommonDispersion <- estimateCommonDisp(NormFactors, verbose=TRUE)");
        cmdSet.add("TagwiseDispersion <- estimateTagwiseDisp(CommonDispersion, trend=\"none\")");
        cmdSet.add("");
        cmdSet.add("ExactTestTagDisp <- exactTest(TagwiseDispersion)");
        cmdSet.add("");
        cmdSet.add("tTags <- topTags(ExactTestTagDisp, n=Inf)");
        cmdSet.add("");
        cmdSet.add("");
        cmdSet.add("resultFile<-\"" + deResultsFile + "\"");
        cmdSet.add("write.table(tTags[tTags$table$PValue<=" + this.getpValue() + ", ], file=\"" + deResultsFile + "\", " + " sep=\",\", row.names=TRUE)");
        
        cmdSet.add("deTags<-rownames(tTags[tTags$table$PValue<=" + this.getpValue()+ ", ])");
        cmdSet.add("");
        cmdSet.add("write.table(cpm(TagwiseDispersion)[deTags,], file=\"" + deCountsBySampleFile + "\", " + " sep=\",\", row.names=TRUE)");
        cmdSet.add("");
        cmdSet.add("write.table(summary(decideTestsDGE(ExactTestTagDisp, p=" + this.getpValue() + ", adjust=\"BH\")), file=\"" + deSummaryFile + "\", " + " sep=\",\")");
        cmdSet.add("");
        
        
        dePlotBCVfile           = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + DE_BCV_PLOT_EXTENSION + "png";
        cmdSet.add("png(\"" + dePlotBCVfile + "\", width=" + PLOT_WIDTH + ", height=" + PLOT_HEIGHT + ", units=\"" + PLOT_UNITS + "\", res=" + PLOT_RESOLUTION + ")");
        cmdSet.add("plotBCV(TagwiseDispersion)");
        cmdSet.add("dev.off()");
        cmdSet.add("");
        
        dePlotMDSfile           = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + DE_MDS_PLOT_EXTENSION + "png";
        cmdSet.add("png(\"" + dePlotMDSfile + "\", width=" + PLOT_WIDTH + ", height=" + PLOT_HEIGHT + ", units=\"" + PLOT_UNITS + "\", res=" + PLOT_RESOLUTION + ")");
        cmdSet.add("plotMDS(TagwiseDispersion)");
        cmdSet.add("dev.off()");
        cmdSet.add("");

        dePlotSmearfile         = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + DE_SMEAR_PLOT_EXTENSION + "png";
        cmdSet.add("png(\"" + dePlotSmearfile + "\", width=" + PLOT_WIDTH + ", height=" + PLOT_HEIGHT + ", units=\"" + PLOT_UNITS + "\", res=" + PLOT_RESOLUTION + ")");
        cmdSet.add("plotSmear(TagwiseDispersion)");
        cmdSet.add("dev.off()");        
        cmdSet.add("");
        
        try{
            BufferedWriter bwRC = new BufferedWriter(new FileWriter(new File(rScriptFilename)));
                for(String cmd: cmdSet){
                    bwRC.write(cmd + "\n");
                }
            bwRC.close();
        }
        catch(IOException exIO){
            logger.error("error writing generated RScript to file <" + rScriptFilename + ">");
            logger.error(exIO);
            throw new IOException(STEP_ID_STRING + ": error writing generated RScript to file <" + rScriptFilename + ">");
        }
    }
    
    
    /**
     * 
     * @throws IOException 
     */
    private void executeRScript() throws IOException{
        

        ArrayList<String> cmd = new ArrayList<>();

        cmd.add("/usr/local/bin/Rscript");
        cmd.add(rScriptFilename);


        String cmdRunRScript = StringUtils.join(cmd, " ");
        cmdRunRScript = cmdRunRScript.replace(FILESEPARATOR + FILESEPARATOR, FILESEPARATOR);
        logger.info("Rscript command:\t" + cmdRunRScript);

        try{
            Runtime rt = Runtime.getRuntime();
            Process proc = rt.exec(cmdRunRScript);
            BufferedReader brAStdin  = new BufferedReader(new InputStreamReader(proc.getInputStream()));
            BufferedReader brAStdErr = new BufferedReader(new InputStreamReader(proc.getErrorStream()));

                String line = null;
                logger.info("<OUTPUT>");
                while ((line = brAStdin.readLine()) != null)
                    logger.info(line);
                logger.info("</OUTPUT>");

                logger.info("<ERROR>");
                while ( (line = brAStdErr.readLine()) != null){
                    logger.info(line);
                }
                // need to parse the output from Bowtie to get the mapping summary
                logger.info("</ERROR>");

                int exitVal = proc.waitFor();            
                logger.info("Process exitValue: " + exitVal);   

            brAStdin.close();
            brAStdErr.close();        
        }
        catch(IOException exIO){
            logger.info("error executing RScript command\n" + cmdRunRScript);
            logger.info(exIO);
            throw new IOException(STEP_ID_STRING + ": error executing RScript command\n" + cmdRunRScript);
            
        }
        catch(InterruptedException exIE){
            logger.info("error executing RScript command\n" + exIE);            
            logger.info(exIE);
            throw new IOException(STEP_ID_STRING + ": error executing RScript command\n" + cmdRunRScript);
        }
    }
    
    
    /**
     * Verify Input Data 
     * 
     * @throws IOException
     */        
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        if (new File(gffFileMirBase).exists()==false){
            logger.error("no annotation file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + gffFileMirBase + ">");
            throw new IOException("no annotation file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + gffFileMirBase + ">");
        }
                
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        if (new File(gffFileMirBase).exists()==false){
            logger.error("no fasta file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + faFileMirBase + ">");
            throw new IOException("no fasta file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + faFileMirBase + ">");
        }
                
//        Iterator itSD = this.getStepInputData().getSampleData().iterator();
//        while (itSD.hasNext()){
//            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
//            String  miRCountsFile  = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MIR_COUNTS_EXTENSION);            
//            if ((new File(miRCountsFile)).exists()==false){
//                throw new IOException(STEP_ID_STRING + ": counts file <" + sampleData + "> does not exist");
//            }
//            
//        }
    }
    
    
    
    
    /**
     * generate sample configuration data so the user can see what can be
     * specified
     *
     * @return
     */
    @Override
    public HashMap generateExampleConfigurationData() {

        logger.info(STEP_ID_STRING + ": generate example configuration data");

        HashMap configData = new HashMap();

        configData.put(ID_REF_GENOME, "hsa");
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_PVALUE, 0.05);

        return configData;
    }


    
    
    
    @Override
    public void verifyOutputData(){
        
    }

    
    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }
    
    
    
    /**
     * @return the pValue
     */
    public double getpValue() {
        return pValue;
    }

    /**
     * @param pValue the pValue to set
     */
    public void setpValue(double pValue) {
        this.pValue = pValue;
    }

    /**
     * @return the miRBaseRelease
     */
    public String getMiRBaseRelease() {
        return miRBaseRelease;
    }

    /**
     * @param miRBaseRelease the miRBaseRelease to set
     */
    public void setMiRBaseRelease(String miRBaseRelease) {
        this.miRBaseRelease = miRBaseRelease;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return ReferenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.ReferenceGenome = ReferenceGenome;
    }


}
