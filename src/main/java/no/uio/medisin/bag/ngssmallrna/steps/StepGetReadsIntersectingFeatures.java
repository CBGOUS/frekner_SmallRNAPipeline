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
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.core.mirna.MiRNA;
import no.uio.medisin.bag.core.mirna.MiRNASet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;

import java.nio.file.Paths;
import no.uio.medisin.bag.core.mapping.MappedRead;
import no.uio.medisin.bag.core.mapping.SAMEntry;
import no.uio.medisin.bag.core.sequence.Strand;
import org.apache.commons.io.FilenameUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   parse SAM file to extract and process the reads that overlap the specified
 *   feature list
 * 
 *   Input is a SAM file of reads and a BED/GFF3 file of feature locations
 *   the SAM files were generated from mapping collapsed FastA reads, where
 *   the header contains the count information
 *   for example, 
 *      >2-4556
 *   means this was the second most abundant read, and it occurred 4556 times
 *   See StepCollapseReads for more details
 * 
 *   Output is 
 *      (i) a tab delimited file of feature counts and 
 *      (ii) a tab delimited file of reads that didn't perfectly align with the features
 * 
 * @author sr
 */

public class StepGetReadsIntersectingFeatures extends NGSStep implements NGSBase{


    
    static Logger                   logger                          = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;
    
    public  static final String     STEP_ID_STRING                  = "GetReadsIntersectingFeatures";
    private static final String     ID_BLEED                        = "bleed";
    private static final String     ID_ISOMIRS                      = "analyzeIsomirs";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_BASELINE                     = "baselinePercent";
    
    private static final String     PARAM_SAMFILE_SHORT             = "-S";
    private static final String     PARAM_SAMFILE_LONG              = "--sam_file";
    private static final String     PARAM_BEDFILE_SHORT             = "-B";
    private static final String     PARAM_BEDFILE_LONG              = "--bed_file";
    private static final String     PARAM_GFFFILE_SHORT             = "-G";
    private static final String     PARAM_GFFFILE_LONG              = "--gff_file";
    private static final String     PARAM_BLEED_SHORT               = "-b";
    private static final String     PARAM_BLEED_LONG                = "--bleed_file";
    private static final String     PARAM_FEATURELIST_SHORT         = "-F";
    private static final String     PARAM_FEATURELIST_LONG          = "--feature_list";
        

    private static final String     INFILE_EXTENSION                = ".trim.clp.gen.sam";
    private static final String     RAW_INPUT_EXTENSION             = ".fastq.gz";
    private static final String     MIRCOUNTS_EXTENSION             = ".trim.clp.gen.mircounts.tsv";
    

    private List<MiRNA>             miRNAHitList;
    MiRNASet                        mirBaseSet                      = new MiRNASet();           
    private ArrayList<MappedRead>   mappedReads                     = new ArrayList<>();
    
    
    private int                     locationBleed                   = 2;
    private String                  referenceGenome                 = "";
    private int                     baselinePercent                 = 5;
    private String                  featureReferenceFile            = "";
    private String []               featureList                     = new String[] {"ALL"};
    private String                  samFileListString               = "";
    private String                  featureListString               = "";
            
    



    public StepGetReadsIntersectingFeatures(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepGetReadsIntersectingFeatures(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    


    @Override
    public String shortStepDescription(){
      return "Find Reads in SAM file Intersecting specified Features ";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "+" + "-".repeat(3) + "+" +
              "parse SAM file to find the reads that intersect the specified features. "
              + "This is primarily intended for counting reads mapping to miRNAs, "
              + "so a 'bleed' variable can be specified to count reads that extend"
              + "past the 5' and 3' ends of the feature." + System.lineSeparator().repeat(2) 
              + "Required Parameters:" + System.lineSeparator().repeat(2)
              + "to run this step, you need specify:" + System.lineSeparator()
              + "either:" + System.lineSeparator()
              + " a feature file in BED format (" + PARAM_BEDFILE_SHORT + "/" + PARAM_BEDFILE_LONG  + ")" + System.lineSeparator()
              + "or:" + System.lineSeparator()
              + " a feature file in GFF format (" + PARAM_GFFFILE_SHORT + "/" + PARAM_GFFFILE_LONG + ")" + System.lineSeparator()
              
              + "Optional Parameters:" + System.lineSeparator().repeat(2)
              + " a list of SAM files (" + PARAM_SAMFILE_SHORT + "/" + PARAM_SAMFILE_LONG+ ")" + System.lineSeparator()              
              + " This is to allow you to specify a subset of SAM files in the InputFolder" + System.lineSeparator()
              + " If you don't specify this parameter, all SAM files specified in the Sample File" + System.lineSeparator()
              + " will be processed" + System.lineSeparator().repeat(2)
              + " a bleed parameter (" + PARAM_BLEED_SHORT + "/" + PARAM_BLEED_LONG+ ")" + System.lineSeparator()              
              + " This is to allow you to count features such as isomiRs, which" + System.lineSeparator()
              + " extend beyond the specified start and stop position." + System.lineSeparator()
              + " For example, setting this to 3 will count reads that are within +/3 nt " + System.lineSeparator()
              + " of the start and stop position of a feature" + System.lineSeparator().repeat(2)
              + " a list of features to check (" + PARAM_FEATURELIST_SHORT + "/" + PARAM_FEATURELIST_LONG + ")" + System.lineSeparator()              
              + " for example, setting this to 'miRNA' will ignore 'miRNA_primary_transcript' entries " + System.lineSeparator()
              + " in a miRBase GFF annotation file when counting reads. " + System.lineSeparator()
              + " Setting to 'ALL' will check all features. " + System.lineSeparator()
              + "Example parameter entry:" + System.lineSeparator() +
              "-b=3, -G=hsa22.1.gff3, -S=sam1.SAM sam2.SAM" + System.lineSeparator().repeat(3);
    }

    @Override
    public void parseStepParameters() throws Exception{
      /*
      here we are looking for BED or GFF file
      and optionally
      Bleed and SAMFile List
      e.g. "-b 3, -G -S"
      */
      logger.info("-".repeat(LOGINDENT) + STEP_ID_STRING + ": parseStepParameters");
      if (this.getStepInputData().getParameterString()==null)
          return;
        
      Boolean bedFile = false;
      Boolean gffFile = false;
      String params[] = this.getStepInputData().getParameterString().split(",");
      for (String param : params){
            String key = param.split("=")[0].trim();
            String value = param.split("=")[1].trim();
            if(key.equals(PARAM_BEDFILE_SHORT) || key.equals(PARAM_BEDFILE_LONG)){
              bedFile = true;     
              this.featureReferenceFile = value.trim();
              logger.info("-".repeat(LOGINDENT) + "--feature reference file set to <" + this.featureReferenceFile + "> ");
            }else if(key.equals(PARAM_GFFFILE_SHORT) || key.equals(PARAM_GFFFILE_LONG)){
              gffFile = true;     
              this.featureReferenceFile = value.trim();
              logger.info("-".repeat(LOGINDENT) + "--feature reference file set to <" + this.featureReferenceFile + "> ");
            }else if(key.equals(PARAM_FEATURELIST_SHORT) || key.equals(PARAM_FEATURELIST_LONG)){ 
              this.featureListString = value.trim();
              this.featureList = featureListString.split(" ");
              logger.info("-".repeat(LOGINDENT) + "--feature list is <" + this.featureListString + "> ");
            }else if(key.equals(PARAM_BLEED_SHORT) || key.equals(PARAM_BLEED_LONG)){
              try{
                String chk = checkParameter("Integer", ID_BLEED, Integer.toString(this.locationBleed), "0", "NA", logger);
                if(chk!=null)
                   this.locationBleed = Integer.parseInt(value.trim());
                   logger.info("-".repeat(LOGINDENT) + "--location bleed set to <" + this.locationBleed + "> ");
              }catch(Exception exEx){
                logger.error("bleed value is out of range <" + this.locationBleed + ">. It must be an integer > 0");         
                logger.error(exEx);
                throw new IllegalArgumentException("bleed value is out of range <" + this.locationBleed + ">. It must be an integer > 0");         
              }
              
            }else if(key.equals(PARAM_SAMFILE_SHORT) || key.equals(PARAM_SAMFILE_LONG)){
              samFileListString = value.trim();
              // here, need to replace the existing file list specified in the run sample file
              // the parameter string will be in the format 'file.sam file2.sam'
              this.getStepInputData().getSampleData().clear();
              String samFiles[] = samFileListString.split(" ");
              for (String samFile: samFiles){
                this.getStepInputData().getSampleData().add(new SampleDataEntry(samFile, "", "", "", "", ""));
              }
            }
      }
      if(bedFile && gffFile){
        logger.error("-".repeat(LOGINDENT) + "Both BED file and GFF files were specified. You can only specify one");
        throw new RuntimeException("Both BED file and GFF files were specified. You can only specify one");
      }
      
    }
    
    
    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * if HashMap is null, then the data entry is probably missing from the YAML file
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception{
        logger.info("-".repeat(LOGINDENT) + STEP_ID_STRING + ": verify configuration data");
        if(configData == null){
          logger.warn("-".repeat(LOGINDENT) + "--: Configuration data is null." + System.lineSeparator()
                  + "-".repeat(LOGINDENT) + "--:Entry is probably missing from the configuration file or the Tag <" + STEP_ID_STRING + "> is misspelled");          
          return;          
        }
          
        if(configData.get(ID_BLEED)==null) {
            logger.error("-".repeat(LOGINDENT) + "--<" + ID_BLEED + "> : Missing Definition in Configuration File");
            throw new NullPointerException("--<" + ID_BLEED + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_ISOMIRS)==null) {
            logger.error("-".repeat(LOGINDENT) + "--<" + ID_ISOMIRS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_ISOMIRS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("-".repeat(LOGINDENT) + "--<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("-".repeat(LOGINDENT) + "--<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_BASELINE)==null) {
            logger.error("-".repeat(LOGINDENT) + "--<" + ID_BASELINE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
        }
        

        String chk;

        chk = checkParameter("Integer", ID_BASELINE, Integer.toString((Integer)configData.get(ID_BASELINE)), "0", "100", logger);
        if(chk!=null)
            this.setBaselinePercent((Integer) configData.get(ID_BASELINE));
        
        chk = checkParameter("Integer", ID_BLEED, Integer.toString((Integer)configData.get(ID_BLEED)), "0", "NA", logger);
        if(chk!=null)
            this.setLocationBleed((Integer) configData.get(ID_BLEED));
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error("-".repeat(LOGINDENT) + ID_REF_GENOME + "--<" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

       
        

        logger.info("-".repeat(LOGINDENT) + "passed");
    }
    
    
    
    /**
     * count up reads that overlap features specified in the GFF file
     * 
     * 
     * @throws IOException 
     */
    @Override
    public void execute()  throws IOException{
        
        logger.info(System.lineSeparator().repeat(2) + "-".repeat(LOGINDENT) + STEP_ID_STRING + ": execute");                
        
        // this needs to load a feature list, 
        String faFileMirBase = featureReferenceFile.replace("gff3", "mature.fa");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(), featureReferenceFile, faFileMirBase);
        
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("-".repeat(LOGINDENT) + "--created output folder <" + outFolder + "> for results" );
        String samLine = null;
        String samInputFilePath = "";
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        int bleed = this.getLocationBleed();
        int i=0;
        while (itSD.hasNext()){
            i=i+1;
            if (i%100000==0){
                logger.info("i");
            }
            try{
                
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                samInputFilePath = Paths.get(inFolder, sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION)).toString();
                logger.info("-".repeat(LOGINDENT) + "-- processing SAM file <" + samInputFilePath + ">");
                if ((new File(samInputFilePath)).exists()==false){
                  logger.warn("-".repeat(LOGINDENT) + "--SAM file <" + samInputFilePath + "> does not exist, skipping");
                  continue;
                }
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalReadCount = 0;
                int totalKnownMiRNACount=0;
                int totalUnknownMiRNACount=0;
                samLine = null;
                BufferedWriter bwTmp = new BufferedWriter(new FileWriter(new File(samInputFilePath + ".txt")));
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFilePath)));
                    //isomiRList = new ArrayList<>();
                    miRNAHitList = new ArrayList<>();
                    String unknowmiRNAstring="unknowmiRNA\n";
                    while((samLine=brSAM.readLine())!= null){

                        SAMEntry samEntry = new SAMEntry(samLine);
                        if (samEntry.isHeaderLine()) 
                            continue;
                        totalReadCount += Integer.parseInt(samEntry.getqName().split("-")[1]);
                        if (samEntry.isMappedRead()) {
                            mappedReads.add(new MappedRead(samEntry.getStartPos(), samEntry.getEndPos(),
                                    samEntry.getrName(), samEntry.getStrand().toString(), Integer.parseInt(samEntry.getqName().split("-")[1])));
                            if(samEntry.getStrand()==Strand.MINUS){
                              preMatchCount3++;
                            }else{
                              preMatchCount5++;
                            }
                            MiRNA miRNAFeature = this.doesReadOverlapKnownMiRNA(samEntry.getStartPos(), 
                                    samEntry.getEndPos(), samEntry.getrName(), samEntry.getStrand().toString(), bleed);
                            if (miRNAFeature != null){
                                int readCount = Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                                
                                totalKnownMiRNACount += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                                bwTmp.write(samLine + System.lineSeparator());
                                MiRNA miRNAHit = new MiRNA(miRNAFeature);
                                //miRNAHit.setReadCount(Integer.parseInt(samLine.split("\t")[0].split("-")[1]));
                                logger.debug("-".repeat(LOGINDENT) + "--miRNA Hit is <" + miRNAHit.getName() + ">");
                                String name = samLine.split("\t")[0];
                                //String sequence = samLine.split("\t")[9];
                                
                                if(miRNAHitList.contains(miRNAHit)){ 
                                    miRNAHitList.get(
                                            miRNAHitList.indexOf(miRNAHit)).addIsomiR(
                                                    samEntry.getqName(), 
                                                    samEntry.getStartPos(), 
                                                    samEntry.getCigar(), 
                                                    samEntry.getMDString(), 
                                                    samEntry.getSeq(),
                                                    readCount);
                                }
                                else{
                                    miRNAHit.addIsomiR(samEntry.getqName(), 
                                            samEntry.getStartPos(), 
                                            samEntry.getCigar(), 
                                            samEntry.getMDString(), 
                                            samEntry.getSeq(),
                                            readCount);
                                    miRNAHitList.add(miRNAHit);
                                }
                                    
                                if(samEntry.getStrand()==Strand.MINUS){
                                  matchCount3++;
                                }else{
                                  matchCount5++;
                                }

                                List<String> mutations = new ArrayList<>();
                                Matcher match = Pattern.compile("[0-9]+|[a-z]+|[A-Z]+").matcher(samEntry.getMDString().split(":")[2]);
                                String outputString = samLine.split("\t")[0] + ":" + samEntry.getStartPos() + ":" + samEntry.getEndPos() + ":[" + samEntry.getrName() + "]:" + samEntry.getMDString() + ": ";
                                while (match.find()) {
                                    mutations.add(match.group());
                                    outputString = outputString.concat(match.group() + "|");
                                }
                            
                            }
                        }
                    }
                    bwTmp.close();
                    logger.info("-".repeat(LOGINDENT) + "--  total mapped counts = " + totalReadCount);                 
                    /*
                        the following is rather approximate.
                        apparently, for 5,000,000 reads, the lowest detectable by qPCR is 50. so, we divide total counts by 100000
                        there has to be a better way....
                    */
                    Double minCountsForSingleFeature = (double) totalReadCount /100000.0; // <= this is rather approximate
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                    
                    
                    logger.info("-".repeat(LOGINDENT) + "--  write miRNA counts");
                    String basename = FilenameUtils.removeExtension(sampleData.getFastqFile1());
                    
                    String  miRCountsFile  = Paths.get(outFolder, basename + MIRCOUNTS_EXTENSION).toString();
                    
                    BufferedWriter brCounts  = new BufferedWriter(new FileWriter(new File(miRCountsFile)));
                        for(MiRNA miR: this.mirBaseSet.getMiRBaseMiRNAList()){
                            if(miRNAHitList.contains(miR)){
                                for(MiRNA miRHit: this.miRNAHitList){
                                    if(miRHit.equals(miR)){
                                        brCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + miRHit.getTotalCounts() + "\n");
                                        break;
                                    }
                                }
                            }
                            else{
                                brCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + 0 + "\n");                                
                            }
                        }
                    brCounts.close();
 
                brSAM.close();      
                logger.info("-".repeat(LOGINDENT) + "--  completed processing SAM file" + System.lineSeparator().repeat(3));
                
                
            }
            catch(IOException ex){
                logger.error("-".repeat(LOGINDENT) + "--error processing sample <" + samInputFilePath + ">\n" + ex.toString());
                throw new IOException(STEP_ID_STRING + ": error processing sample <" + samInputFilePath + ">");
            }
            catch(ArrayIndexOutOfBoundsException exBnd){
                logger.error("-".repeat(LOGINDENT) + "--error parsing line " + samLine);
                logger.error(exBnd);
                throw new IOException(STEP_ID_STRING + ": error processing sample <" + samInputFilePath + ">: samLine was \n" + samLine);
            }
        }
        
        
        logger.info("-".repeat(LOGINDENT) + "--" + STEP_ID_STRING + ": completed");
 
    }
    
    
    
    
    /**
     * Does the read sufficiently overlap a defined miRNA entry?
     * 
     * @param start
     * @param stop
     * @param chr
     * @param strand
     * @param bleed     int : specifies how much a read can 'miss' an entry
     *                        and still be counted
     * 
     * @return MiRNAFeature
     */
    public MiRNA doesReadOverlapKnownMiRNA(int start, int stop, String chr, String strand, int bleed){
        
        for(MiRNA miRBaseEntry: this.mirBaseSet.getMiRBaseMiRNAList()){
            if (miRBaseEntry.chromosomeMatch(chr)){
                if(strand.equals(miRBaseEntry.getStrand())){


                    
                    if( java.lang.Math.abs(start - miRBaseEntry.getStartPos()) <= bleed){

                        if( java.lang.Math.abs(stop - miRBaseEntry.getEndPos()) <= bleed){
                            return miRBaseEntry;                            
                        }

                    }
                    
                }
            }
        }
        
        return null;
        
    }
    
    
    
    
    /**
     * Verify Input Data for parsing SAM file for miRNAs
     * 
     */        
    @Override
    public void verifyInputData() throws IOException{
        
        logger.info("verify input data");        
        this.setPaths();
        
        if (new File(this.featureReferenceFile).exists()==false){
            logger.error("The specified Feature Reference File :<" + this.featureReferenceFile + "> not found");
            throw new IOException("The specified Feature Reference File :<" + this.featureReferenceFile + "> not found");
        }
                                
        // check the SAM files exist
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            if (sampleData.getFastqFile1()==null) throw new IOException("no Fastq1 file specified");
            String samInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);            
            
            if ((new File(samInputFile)).exists()==false){
                logger.warn(STEP_ID_STRING + ": SAM file <" + samInputFile + "> does not exist");
                // (This was changed to a warning, rather than throwing an exception, because the user
                // may only want to analyze a subset of results)  23/1/2020
                //throw new IOException(STEP_ID_STRING + ": SAM file <" + samInputFile + "> does not exist");
            }
            if (samInputFile.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + samInputFile + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
        }

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
        configData.put(ID_BLEED, 2);
        configData.put(ID_BASELINE, 5);
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_ISOMIRS, true);

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
     * @return the locationBleed
     */
    public int getLocationBleed() {
        return locationBleed;
    }

    /**
     * @param locationBleed the locationBleed to set
     */
    public void setLocationBleed(int locationBleed) {
        this.locationBleed = locationBleed;
    }

    /**
     * @return the ReferenceGenome
     */
    public String getReferenceGenome() {
        return referenceGenome;
    }

    /**
     * @param ReferenceGenome the ReferenceGenome to set
     */
    public void setReferenceGenome(String ReferenceGenome) {
        this.referenceGenome = ReferenceGenome;
    }

    /**
     * @return the baselinePercent
     */
    public int getBaselinePercent() {
        return baselinePercent;
    }

    /**
     * @param baselinePercent the baselinePercent to set
     */
    public void setBaselinePercent(int baselinePercent) {
        this.baselinePercent = baselinePercent;
    }

    /**
     * @return the bedReferenceFile
     */
    public String getBedReferenceFile() {
      return featureReferenceFile;
    }

    /**
     * @param bedReferenceFile the bedReferenceFile to set
     */
    public void setBedReferenceFile(String bedReferenceFile) {
      this.featureReferenceFile = bedReferenceFile;
    }

}
