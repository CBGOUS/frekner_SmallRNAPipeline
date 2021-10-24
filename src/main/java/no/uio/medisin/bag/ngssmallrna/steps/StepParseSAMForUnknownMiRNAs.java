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
import no.uio.medisin.bag.core.mirna.IsomiRSet;
import no.uio.medisin.bag.core.mirna.MiRNA;
import no.uio.medisin.bag.core.mirna.MiRNASet;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * search for potential miRNA that are not present in miRBase
 * 
 * @author joey
 */
public class StepParseSAMForUnknownMiRNAs extends NGSStep{
    
    static Logger                   logger                          = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;
    
    public  static final String     STEP_ID_STRING                  = "ParseSAMForUnknownMiRNAs";
    private static final String     ID_BLEED                        = "bleed";
    private static final String     ID_ISOMIRS                      = "analyzeIsomirs";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";
    private static final String     ID_HOST                         = "host";
    private static final String     ID_BASELINE                     = "baselinePercent";
        

    private static final String     INFILE_EXTENSION                = ".trim.clp.gen.sam";
    private static final String     RAW_INPUT_EXTENSION             = ".fastq.gz";
    private static final String     UNKNOWN_SMALL_RNA_EXTENSION     = ".fastq.gz";
//    private static final String     ISOMIR_SUMMARY_EXTENSION        = ".trim.clp.gen.iso_summary.tsv";
//    private static final String     ISOMIR_PRETTY_EXTENSION         = ".trim.clp.gen.iso_pretty.tsv";
//    private static final String     MIRCOUNTS_EXTENSION             = ".trim.clp.gen.mircounts.tsv";
    

    private List<MiRNA>             miRNAHitList;
    private ArrayList<IsomiRSet>    isomiRList;
    MiRNASet                        mirBaseSet                      = new MiRNASet();           
    
    private List<MiRNA>             unknownMiRNAsList;
    private ArrayList<IsomiRSet>    unknownIsomiRList;
    
    private int                     locationBleed                   = 2;
    private Boolean                 analyzeIsomirs                  = false;
    private int                     miRBaseRelease                  = 20;
    private String                  referenceGenome                 = "";
    private int                     baselinePercent                 = 5;
    
    
    public StepParseSAMForUnknownMiRNAs(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepParseSAMForUnknownMiRNAs(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    

    @Override
    public String shortStepDescription(){
      return "search for potential miRNA that are not present in miRBase";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "search for potential miRNA that are not present in miRBase.\n"
              + "Input is a SAM file.\n";
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    
    /**
     * Verify Input Data for parsing SAM file for unknown miRNAs
     * 
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
                
        // check the SAM files exist
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        while (itSD.hasNext()){
            SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
            if (sampleData.getFastqFile1()==null) throw new IOException("no Fastq1 file specified");
            String samInputFile = inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION);            
            
            if ((new File(samInputFile)).exists()==false){
                throw new IOException(STEP_ID_STRING + ": SAM file <" + samInputFile + "> does not exist");
            }
            if (samInputFile.toUpperCase().endsWith(INFILE_EXTENSION.toUpperCase())==false)
            {
                throw new IOException(STEP_ID_STRING + ": incorrect file extension for input file <" 
                  + samInputFile + ">.  \n" 
                  + "should have <" + INFILE_EXTENSION + "> as extension");
            }
            
        }

    }

    @Override
    public void verifyOutputData() throws IOException {
        
    }

    /**
     * in this method we are simply checking that the configuration file 
     * has all the entries we need. We dont check if the values are acceptable
     * that is the role of the NGSStep.
     * 
     * @param configData
     * @throws Exception 
     */
    @Override
    public void parseConfigurationData(HashMap configData) throws Exception {
        logger.info(STEP_ID_STRING + ": verify configuration data");
        
        if(configData.get(ID_BLEED)==null) {
            logger.error("<" + ID_BLEED + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BLEED + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_ISOMIRS)==null) {
            logger.error("<" + ID_ISOMIRS + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_ISOMIRS + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_MIRBASE_VERSION)==null) {
            logger.error("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_MIRBASE_VERSION + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_HOST)==null) {
            logger.error("<" + ID_HOST + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_HOST + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_BASELINE)==null) {
            logger.error("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
        }
        

        String chk;
        chk = checkParameter("Integer", ID_MIRBASE_VERSION, Integer.toString((Integer)configData.get(ID_MIRBASE_VERSION)), "0", "NA", logger);
        if(chk!=null)
            this.setMiRBaseRelease((Integer) configData.get(ID_MIRBASE_VERSION));

        chk = checkParameter("Integer", ID_BASELINE, Integer.toString((Integer)configData.get(ID_BASELINE)), "0", "100", logger);
        if(chk!=null)
            this.setBaselinePercent((Integer) configData.get(ID_BASELINE));
        
        chk = checkParameter("Integer", ID_BLEED, Integer.toString((Integer)configData.get(ID_BLEED)), "0", "NA", logger);
        if(chk!=null)
            this.setLocationBleed((Integer) configData.get(ID_BLEED));
        

        this.setReferenceGenome((String) configData.get(ID_HOST));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_HOST + " <" + configData.get(ID_HOST) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_HOST + " <" + configData.get(ID_HOST) + "> must be a 3 letter string");            
        }

        chk = checkParameter("Boolean", ID_ISOMIRS, Boolean.toString((Boolean)configData.get(ID_ISOMIRS)), "0", "NA", logger);
        if(chk!=null)
            this.setAnalyzeIsomirs((Boolean) configData.get(ID_ISOMIRS));        

        logger.info("passed");
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

        configData.put(ID_HOST, "hsa");
        configData.put(ID_BLEED, 2);
        configData.put(ID_BASELINE, 5);
        configData.put(ID_MIRBASE_VERSION, 20);
        configData.put(ID_ISOMIRS, true);

        return configData;   
    }

    /**
     * count up reads that overlap features specified in the GFF file
     * 
     * @throws IOException 
     */
    @Override
    public void execute() throws IOException {
        logger.info(STEP_ID_STRING + ": execute");                
        
    
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        String faFileMirBase = gffFileMirBase.replace("gff3", "mature.fa");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);
        
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        String samLine = null;
        String samInputFile = "";
//        String countSummary="#name,totalReads,knownMiRNAcount,unknownMiRNAcount \n";
        String countSummary="";
        Iterator itSD = this.getStepInputData().getSampleData().iterator();
        int bleed = this.getLocationBleed();
        while (itSD.hasNext()){
            try{
                
                SampleDataEntry sampleData = (SampleDataEntry)itSD.next();
                
                samInputFile = this.cleanPath(inFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                logger.info(sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, INFILE_EXTENSION));
                int matchCount5 = 0;
                int matchCount3 = 0;
                int preMatchCount5 = 0;
                int preMatchCount3 = 0;
                int totalCounts = 0;
                int totalKnownMiRNACount=0;
                int totalUnknownMiRNACount=0;
                samLine = null;
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));
                    isomiRList = new ArrayList<>();
                    miRNAHitList = new ArrayList<>();
                    unknownMiRNAsList = new ArrayList<>();
                    unknownIsomiRList = new ArrayList<>();
                    String unknowmiRNAstring="";
                    while((samLine=brSAM.readLine())!= null){
                        /*
                            1   QNAME	   Query template NAME
                            2   FLAG	   bitwise FLAG
                            3   RNAME	   Reference sequence NAME
                            4   POS	   1-based leftmost mapping POSition
                            5   MAPQ	   MAPping Quality
                            6   CIGAR	   CIGAR string
                            7   RNEXT	   Ref. name of the mate/next read
                            8   PNEXT	   Position of the mate/next read
                            9   TLEN	   observed Template LENgth
                            10  SEQ	   segment SEQuence
                            11  QUAL	   ASCII of Phred-scaled base QUALity+33
                        
                        */
                        if(samLine.startsWith("@")) continue;
                        
                        totalCounts += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                        if(samLine.split("\t")[1].equals("16") || samLine.split("\t")[1].equals("0")){
                            String strand = "";
                            if (samLine.split("\t")[1].equals("16")) {
                                strand = "-";
                                preMatchCount3++;
                            }
                            else{
                                strand = "+";
                                preMatchCount5++;
                            }
                            
                            int startPos = Integer.parseInt(samLine.split("\t")[3]);
                            String cigarStr = samLine.split("\t")[5].replace("M", "").trim();
                            int endPos = startPos + Integer.parseInt(cigarStr);
                            String chr = samLine.split("\t")[2].trim();
                            String mdString = samLine.split("\t")[12];
                            String recordName = samLine.split("\t")[0];
                            String recordSeq= samLine.split("\t")[9];
 
                            MiRNA miRNAFeature = this.doesReadOverlapKnownMiRNA(startPos, endPos, chr, strand, bleed);
                            if (miRNAFeature != null){
                                totalKnownMiRNACount += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
//                                MiRNA miRNAHit = new MiRNA(miRNAFeature);
                                //logger.debug(miRNAHit.getName());
//                                String name = samLine.split("\t")[0];
//                                String sequence = samLine.split("\t")[9];
  
                            }
                            else {
                                logger.debug("unknowmiRNA --- " + recordName+"---"+chr +"---"+ startPos+"---"+ endPos+"---"+ strand+"---"+ recordSeq);
                                if (Integer.parseInt(recordName.split("-")[1]) > 10){
                                    if (unknowmiRNAstring.length() > 2){
                                        unknowmiRNAstring=unknowmiRNAstring.concat(recordName+","+chr +","+cigarStr+","+ startPos+","+ endPos+","+ strand+","+ recordSeq+"\n");
                                    }
                                    else{
                                        unknowmiRNAstring=(recordName+","+chr +","+cigarStr+","+ startPos+","+ endPos+","+ strand+","+ recordSeq+"\n");
                                    }
                                MiRNA unknowMiRNA = new MiRNA(recordName, chr, startPos, endPos, strand, "unknown", "unknown", recordSeq);
                                logger.debug("locating error --> " + recordName);
//                                checkOverlapUnknowMiRNA(unknowMiRNA,cigarStr,mdString);
                                }

                                totalUnknownMiRNACount += Integer.parseInt(samLine.split("\t")[0].split("-")[1]);

                            }                        
                        }
                    }
                    logger.info("  total mapped counts = " + totalCounts);
                    logger.debug("totalCount"+ totalCounts+" --- "+ totalKnownMiRNACount +" -- "+ totalUnknownMiRNACount);
                    countSummary+=(sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "")+","+totalCounts+","+ totalKnownMiRNACount +","+ totalUnknownMiRNACount+"\n");
                    /*
                        the following is rather approximate.
                        apparently, for 5,000,000 reads, the lowest detectable by qPCR is 50. so, we divide total counts by 100000
                        there has to be a better way....
                    */
                    Double minCountsForSingleFeature = (double) totalCounts /100000.0; // <= this is rather approximate
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                    
                    if(analyzeIsomirs){
                        logger.info("  calculate isomiR dispersions");
                        for(MiRNA miRHit: miRNAHitList){
                            if (miRHit.getTotalCounts() > minCountsForSingleFeature.intValue()){
                                ArrayList isomirPtsAsHash = miRHit.characterizeIsomiRs(this.getBaselinePercent());
                                this.isomiRList.add(new IsomiRSet(miRHit.getMimatID(), sampleData.getNote(), sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ""), isomirPtsAsHash));
                            }
                        }




                    //analysis unknown miRNA isoform
                        logger.debug("  write isomiRs for unknow miRNAs");

//                        String  isoUnknownDetailsFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".unknown-iso_summary.tsv");
//                        String  isoUnknownPrettyFile  = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".unknown-iso_pretty.tsv");

//                        BufferedWriter brDetailsUnknown = new BufferedWriter(new FileWriter(new File(isoUnknownDetailsFile)));
//                        BufferedWriter brPrettyUnknown  = new BufferedWriter(new FileWriter(new File(isoUnknownPrettyFile)));
//                        logger.debug("miRUnknow  "+unknownMiRNAsList.size());
//                            for(MiRNA miRUnknow: this.unknownMiRNAsList){
////                                if (miRUnknow.getTotalCounts() > minCountsForSingleFeature.intValue()){
//                                    logger.debug("miRUnknow  "+miRUnknow.getName());
//                                    brDetailsUnknown.write(miRUnknow.reportIsomiRs(this.getBaselinePercent(), minCountsForSingleFeature.intValue()));
//                                    brPrettyUnknown.write(miRUnknow.prettyReportIsomiRs(this.getBaselinePercent(), minCountsForSingleFeature.intValue()));
////                                }
//                            }
//                        brDetailsUnknown.close();
//                        brPrettyUnknown.close();


//                        String subfolder = outFolder+"/unknownMiRNAs";
////                        Boolean ufA = new File(subfolder).mkdir();       
////                        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
                        String  isoUnknownDetailsFilename = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".unknowntable.tsv");
                        BufferedWriter brtableUnknown  = new BufferedWriter(new FileWriter(new File(isoUnknownDetailsFilename)));        
                        brtableUnknown.write(unknowmiRNAstring);
                        brtableUnknown.close();
                        
                    
                    //end of analysis unknown miRNA isoform 
                    }

                     // write unknown miRNA counts
//                logger.info("  write unknown miRNA counts");
//
//                String miRUnknownCountsFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ".unknow.mircount.tsv");
//
//                BufferedWriter brUnknownCounts = new BufferedWriter(new FileWriter(new File(miRUnknownCountsFile)));
//                for (MiRNA miR : this.mirBaseSet.getMiRBaseMiRNAList()) {
//                    if (unknownMiRNAsList.contains(miR)) {
//                        for (MiRNA miRHit : this.unknownMiRNAsList) {
//                            if (miRHit.equals(miR)) {
//                                brUnknownCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + miRHit.getTotalCounts() + "\n");
//                                break;
//                            }
//                        }
//                    } else {
//                        brUnknownCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + 0 + "\n");
//                    }
//                }
//                brUnknownCounts.close();

                // end of  write unknown miRNA counts     
                    
                brSAM.close();      
                logger.info("  completed processing SAM file\n\n");
                
                
            }
            catch(IOException ex){
                logger.error("error processing sample <" + samInputFile + ">\n" + ex.toString());
                throw new IOException(STEP_ID_STRING + ": error processing sample <" + samInputFile + ">");
            }
            catch(ArrayIndexOutOfBoundsException exBnd){
                logger.error("error parsing line " + samLine);
                logger.error(exBnd);
                throw new IOException(STEP_ID_STRING + ": error processing sample <" + samInputFile + ">: samLine was \n" + samLine);
            }
        }
        
        if(analyzeIsomirs){
            String dispersionFile   = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".disp.tsv";
            String summaryFile      = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".disp.summary.tsv";
            logger.info("write dispersions to file <" + dispersionFile + ">");
            try{
                BufferedWriter bwDp = new BufferedWriter(new FileWriter(new File(dispersionFile)));
                BufferedWriter bwSm = new BufferedWriter(new FileWriter(new File(summaryFile)));   
                    bwSm.write(IsomiRSet.printSummaryHeader());
                    for(IsomiRSet isomiRset: isomiRList){
                        isomiRset.calcDistParameters();
                        bwSm.write(isomiRset.printSummary());
                        bwDp.write(isomiRset.tabReportIsomiRSet());                
                    }
                bwSm.close();
                bwDp.close();
            }
            catch(IOException exIO){
                logger.info("error writing isomiR dispersion File <" + dispersionFile + ">\n" + exIO);
                throw new IOException(STEP_ID_STRING + "error writing isomiR dispersion File <" + dispersionFile + ">");
            }
           
        }
        
        // write summary of known and unknown miRNA counts 
        logger.info("  write summary of known and unknown miRNA counts ");
        logger.info(countSummary);
        String countSummaryFile = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".summaryall.mircount.tsv";
        BufferedWriter brsummaryCounts = new BufferedWriter(new FileWriter(new File(countSummaryFile),true));
        brsummaryCounts.write(countSummary);
        brsummaryCounts.close();
        // end of  write summary of known and unknown miRNA counts 

        //raw idea: read read in files, check similar reads inthe same file, idenfity them, name them and go to next
        //change non DATAABLE??
        
        
//        checkumiR();
        
        logger.info(STEP_ID_STRING + ": completed");
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

                    if( java.lang.Math.abs(start - miRBaseEntry.getMiStart()) <= bleed){

                        if( java.lang.Math.abs(stop - miRBaseEntry.getMiEnd()) <= bleed){
                            return miRBaseEntry;                            
                        }

                    }

                }
            }
        }
        
        return null;
        
    }
    
    
    
    
    @Override
    public NGSStepSubclass getStepSubclass(){
        return STEP_SUBCLASS;
    }

    /**
     * @return the miRBaseRelease
     */
    public int getMiRBaseRelease() {
        return miRBaseRelease;
    }

    /**
     * @param miRBaseRelease the miRBaseRelease to set
     */
    public void setMiRBaseRelease(int miRBaseRelease) {
        this.miRBaseRelease = miRBaseRelease;
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
     * @return the analyzeIsomirs
     */
    public Boolean getAnalyzeIsomirs() {
        return analyzeIsomirs;
    }

    /**
     * @param analyzeIsomirs the analyzeIsomirs to set
     */
    public void setAnalyzeIsomirs(Boolean analyzeIsomirs) {
        this.analyzeIsomirs = analyzeIsomirs;
    }
    
    private void checkOverlapUnknowMiRNA(MiRNA tmpUnknowMiRNA, String cigarStr,String mdString ) {
        if (unknownMiRNAsList.size() != 0) {
            logger.debug("locationBleed"+tmpUnknowMiRNA);
            for (MiRNA recordInList : unknownMiRNAsList) {
                logger.debug("recordInList"+recordInList);
                if (java.lang.Math.abs(tmpUnknowMiRNA.getMiStart() - recordInList.getMiStart()) <= locationBleed) {
                    if (java.lang.Math.abs(tmpUnknowMiRNA.getMiEnd() - recordInList.getMiEnd()) <= locationBleed) {
                        unknownMiRNAsList.get(unknownMiRNAsList.indexOf(recordInList)).addIsomiR(tmpUnknowMiRNA.getName(), tmpUnknowMiRNA.getMiStart(), cigarStr, mdString, tmpUnknowMiRNA.getSeq(),0);
                        logger.debug("NO ConcurrentModificationException-1");
                    }
                    else {
//                        tmpUnknowMiRNA.addIsomiR(tmpUnknowMiRNA.getName(), tmpUnknowMiRNA.getStartPos(), cigarStr, mdString, tmpUnknowMiRNA.getSequence());
                        unknownMiRNAsList.add(tmpUnknowMiRNA);
                        logger.debug("NO ConcurrentModificationException-2");
                    }
                } 
                else {
//                    tmpUnknowMiRNA.addIsomiR(tmpUnknowMiRNA.getName(), tmpUnknowMiRNA.getStartPos(), cigarStr, mdString, tmpUnknowMiRNA.getSequence());
                    unknownMiRNAsList.add(tmpUnknowMiRNA);
                    logger.debug("NO ConcurrentModificationException-3");
                }
            }
        } 
        else {
//            tmpUnknowMiRNA.addIsomiR(tmpUnknowMiRNA.getName(), tmpUnknowMiRNA.getStartPos(), cigarStr, mdString, tmpUnknowMiRNA.getSequence());
            unknownMiRNAsList.add(tmpUnknowMiRNA);
            logger.debug("NO ConcurrentModificationException-4");
        }
    }
    
}
