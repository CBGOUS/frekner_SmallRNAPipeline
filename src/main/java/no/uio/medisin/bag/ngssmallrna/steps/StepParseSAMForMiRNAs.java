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
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import no.uio.medisin.bag.core.mirna.IsomiRSet;
import no.uio.medisin.bag.core.mirna.IsomiRString;
import no.uio.medisin.bag.core.mirna.MiRNA;
import no.uio.medisin.bag.core.mirna.MiRNASet;
import no.uio.medisin.bag.core.sequence.SimpleSeq;
import no.uio.medisin.bag.core.sequence.Strand;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import org.apache.commons.lang3.math.NumberUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   parse SAM file to extract and process the miRNA reads to determine isomiR content
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class StepParseSAMForMiRNAs extends NGSStep implements NGSBase{
    
    static Logger                   logger                          = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;
    
    public  static final String     STEP_ID_STRING                  = "ParseSAMForMiRNAs";
    private static final String     ID_BLEED                        = "bleed";
    private static final String     ID_ISOMIRS                      = "analyzeIsomirs";
    private static final String     ID_GROUP_BY_SEED                = "groupBySeedRegion";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_BASELINE                     = "baselinePercent";
        

    private static final String     INFILE_EXTENSION                = ".trim.clp.gen.sam";
    private static final String     RAW_INPUT_EXTENSION             = ".fastq.gz";
    private static final String     ISOMIR_SUMMARY_EXTENSION        = ".trim.clp.gen.iso_summary.tsv";
    private static final String     ISOMIR_PRETTY_EXTENSION         = ".trim.clp.gen.iso_pretty.tsv";
    private static final String     MIRCOUNTS_EXTENSION             = ".trim.clp.gen.mircounts.tsv";
    

    private List<MiRNA>             miRNAHitList;
    private ArrayList<IsomiRSet>    isomiRList;
    MiRNASet                        mirBaseSet                      = new MiRNASet();           
    
    
    private int                     locationBleed                   = 2;
    private Boolean                 analyzeIsomirs                  = false;
    private Boolean                 groupBySeedRegion               = false;        
    private String                  miRBaseRelease                  = "22.1";
    private String                  referenceGenome                 = "";
    private int                     baselinePercent                 = 5;
    



    public StepParseSAMForMiRNAs(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepParseSAMForMiRNAs(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    

    @Override
    public String shortStepDescription(){
      return "parse SAM file to extract and process the miRNA reads to determine isomiR content";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "parse SAM file to extract and process the miRNA reads to determine isomiR content.\n"
              + "TNote: a much richer analysis can be used with the standalone"
              + "Jasmine package.\n"
              + "However, this step is more readily customizable\n";
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
    public void parseConfigurationData(HashMap configData) throws Exception{
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
        if(configData.get(ID_REF_GENOME)==null) {
            logger.error("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_REF_GENOME + "> : Missing Definition in Configuration File");
        }
        if(configData.get(ID_BASELINE)==null) {
            logger.error("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
            throw new NullPointerException("<" + ID_BASELINE + "> : Missing Definition in Configuration File");
        }
        

        String chk;
        /* We shouldn't assume the miRBase version is a number, so remove this check.
        chk = checkParameter("Double", ID_MIRBASE_VERSION, Double.toString((Double)configData.get(ID_MIRBASE_VERSION)), "0", "NA", logger);
        if(chk!=null)
            this.setMiRBaseRelease((Double) configData.get(ID_MIRBASE_VERSION));
        */
        this.setMiRBaseRelease((String)configData.get(ID_MIRBASE_VERSION));
        chk = checkParameter("Integer", ID_BASELINE, Integer.toString((Integer)configData.get(ID_BASELINE)), "0", "100", logger);
        if(chk!=null)
            this.setBaselinePercent((Integer) configData.get(ID_BASELINE));
        
        chk = checkParameter("Integer", ID_BLEED, Integer.toString((Integer)configData.get(ID_BLEED)), "0", "NA", logger);
        if(chk!=null)
            this.setLocationBleed((Integer) configData.get(ID_BLEED));
        

        this.setReferenceGenome((String) configData.get(ID_REF_GENOME));
        if(this.getReferenceGenome().length() !=3 ){
            logger.error(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
            throw new IllegalArgumentException(ID_REF_GENOME + " <" + configData.get(ID_REF_GENOME) + "> must be a 3 letter string");            
        }

        chk = checkParameter("Boolean", ID_ISOMIRS, Boolean.toString((Boolean)configData.get(ID_ISOMIRS)), "0", "NA", logger);
        if(chk!=null)
            this.setAnalyzeIsomirs((Boolean) configData.get(ID_ISOMIRS));        
        
        chk = checkParameter("Boolean", ID_GROUP_BY_SEED, Boolean.toString((Boolean)configData.get(ID_GROUP_BY_SEED)), "0", "NA", logger);
        if(chk!=null)
            this.setGroupBySeedRegion((Boolean) configData.get(ID_GROUP_BY_SEED));        

        logger.info("passed");
    }
    
    
    
    /**
     * count up reads that overlap features specified in the GFF file
     * 
     * @throws IOException 
     */
    @Override
    public void execute()  throws IOException{
        
        logger.info(STEP_ID_STRING + ": execute");                
        
    
        String gffFileMirBase = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        String faFileMirBase = gffFileMirBase.replace(".gff3", "_mature.fa");
        mirBaseSet.loadMiRBaseData(this.getReferenceGenome(), gffFileMirBase, faFileMirBase);
        
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        String samLine = null;
        String samInputFile = "";
        HashMap <String, String> uniqMirNames = new HashMap<>();
        
        Path uniqIsomiRsFilePath = Paths.get(outFolder,this.getStepInputData().getProjectID()+"_uniq_isomirs.tsv");
        logger.info("isomiR information will be written to <" + uniqIsomiRsFilePath.toString() + ">");
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
                    String unknowmiRNAstring="unknowmiRNA\n";
                    int i=0;
                    while((samLine=brSAM.readLine())!= null){
                        i=i+1;
                        if (i%100000==0){
                            logger.info(i);
                        }
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
                            int endPos = startPos + Integer.parseInt(cigarStr)-1;
                            String chr = samLine.split("\t")[2].trim();
                            String mdString = samLine.split("\t")[12];
                            String recordName = samLine.split("\t")[0];
                            String recordSeq= samLine.split("\t")[9];
                            
                            
                            
                            MiRNA miRNAFeature = this.mirBaseSet.doesReadOverlapKnownMiRNA(startPos, endPos, chr, strand, bleed);
                            if (miRNAFeature != null){
                                int readCount = Integer.parseInt(samLine.split("\t")[0].split("-")[1]);
                                totalKnownMiRNACount += readCount;
                                MiRNA miRNAHit = new MiRNA(miRNAFeature);
                                if (miRNAHit.getStrand().equals("-")){
                                  logger.debug("- strand");
                                }
                                //logger.debug(miRNAHit.getName());
                                String readName = miRNAHit.getName();
                                
                                if(startPos!=miRNAHit.getStartPos() || endPos!=miRNAHit.getEndPos()){
                                  readName =  readName + "|iso";
                                }else{
                                  String mdField = mdString.split(":")[2].trim();
                                  if(NumberUtils.isNumber(mdField)){
                                    int ln = Integer.valueOf(mdField);
                                    if(ln!=miRNAHit.getLength()){
                                      readName =readName + "|iso";
                                    }
                                  }else{
                                    readName = readName + "|iso";
                                  }
                                }
    
                                String sequence = samLine.split("\t")[9];
                                
                                if(miRNAHitList.contains(miRNAHit)){           
                                    miRNAHitList.get(miRNAHitList.indexOf(miRNAHit)).addIsomiR(readName, startPos, cigarStr, mdString, sequence, readCount);
                                }
                                else{
                                    miRNAHit.addIsomiR(readName, startPos, cigarStr, mdString, sequence, readCount);
                                    miRNAHitList.add(miRNAHit);
                                }
                                    
                                if(strand.equals("+")) matchCount5++;
                                else matchCount3++;
                                List<String> mutations = new ArrayList<>();
                                Matcher match = Pattern.compile("[0-9]+|[a-z]+|[A-Z]+").matcher(samLine.split("\t")[12].split(":")[2]);
                                String outputString = samLine.split("\t")[0] + ":" + startPos + ":" + endPos + ":[" + chr + "]:" + samLine.split("\t")[12] + ": ";
                                while (match.find()) {
                                    mutations.add(match.group());
                                    outputString = outputString.concat(match.group() + "|");
                                }
                                
                            }
       
                        }
                    }
                    logger.info("  total mapped counts = " + totalCounts);                 
                    /*
                        the following is rather approximate.
                        apparently, for 5,000,000 reads, the lowest detectable by qPCR is 50. so, we divide total counts by 100000.
                        There has to be a better way....
                    */
                    Double minCountsForSingleFeature = (double) totalCounts /100000.0; // <= this is rather approximate
                    logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped");
                    //analyzeIsomirs=false;
                    if(analyzeIsomirs){
                        logger.info("  calculate isomiR dispersions");
                        for(MiRNA miRHit: miRNAHitList){
                            if (miRHit.getTotalIsomiRCounts() > minCountsForSingleFeature.intValue()){
                                ArrayList isomirPtsAsHash = miRHit.characterizeIsomiRs(this.getBaselinePercent());
                                this.isomiRList.add(new IsomiRSet(miRHit.getMimatID(), sampleData.getNote(), sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ""), isomirPtsAsHash));
                            }
                        }



                        logger.info("  write isomiRs");

                        String  isoDetailsFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ISOMIR_SUMMARY_EXTENSION);
                        String  isoPrettyFile  = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, ISOMIR_PRETTY_EXTENSION);

                        BufferedWriter brDetails = new BufferedWriter(new FileWriter(new File(isoDetailsFile)));
                        BufferedWriter brPretty  = new BufferedWriter(new FileWriter(new File(isoPrettyFile)));
                            for(MiRNA miRHit: this.miRNAHitList){
                                if (miRHit.getTotalIsomiRCounts() > minCountsForSingleFeature.intValue()){
                                    logger.debug(miRHit.getName());
                                    brDetails.write(miRHit.reportIsomiRsAsTSVs(this.getBaselinePercent(), minCountsForSingleFeature.intValue()));
                                    brPretty.write(miRHit.prettyReportIsomiRs2(this.getBaselinePercent(), minCountsForSingleFeature.intValue()));
                                }
                            }
                        brPretty.close();
                        brDetails.close();

                    }
                    
                    /*
                    here we have to build a list of unique miR/isomiR names
                    as we will need this information to integrate the reads counts
                    from all the other samples
                    
                    to generate a unique miR name we need the MIMAT id + miR name
                    + startPos/MD String/
                    */
                    logger.info("  write miRNA counts");
                    
                    if(this.getGroupBySeedRegion()){
                      for(MiRNA miRHit: this.miRNAHitList){
                        miRHit.groupIsomiRsBySeedRegion(this.getBaselinePercent());
                      // update the existing isomiR hashmap with the Group By Seed hashmap
                      }

                    }            
                    
                    


                    String  miRCountsFile  = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, MIRCOUNTS_EXTENSION);
                    
                    BufferedWriter brCounts  = new BufferedWriter(new FileWriter(new File(miRCountsFile)));
                        for(MiRNA miR: this.mirBaseSet.getMiRBaseMiRNAList()){
                            if(miRNAHitList.contains(miR)){
                              //miRNAHitList.get(miRNAHitList.indexOf(miRNAHit))
                                for(MiRNA miRHit: this.miRNAHitList){
                                    if(miRHit.equals(miR)){
                                        // write out all the isomiRs and add the names to the global list
                                        for(IsomiRString isoStr: miRHit.isomiRs.values()){
                                          String isoName = "";
                                          if(isoStr.getName().contains("iso")){
                                            // if this is an isomiR, need to write enough info to distinguish it
                                            isoName = isoStr.getName() + "|" + isoStr.getStart() + "|" + isoStr.getMDString();
                                          }else{
                                            isoName = isoStr.getName();
                                          }
                                          brCounts.write(miR.getMimatID() + ":" + isoName + "\t" + isoStr.getCount() + "\n"); 
                                          if(miR.getStrand().equals(Strand.MINUS.toString())){
                                            uniqMirNames.put(miR.getMimatID() + ":" + isoName, SimpleSeq.complement(isoStr.getSeq()));
                                          }
                                          else{
                                            uniqMirNames.put(miR.getMimatID() + ":" + isoName, isoStr.getSeq());
                                          }
                                        }
                                        //brCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + miRHit.getTotalCounts() + "\n");
                                        break;
                                    }
                                }
                            }
                            else{
                                brCounts.write(miR.getMimatID() + ":" + miR.getName() + "\t" + 0 + "\n"); 
                                uniqMirNames.put(miR.getMimatID() + ":" + miR.getName(), SimpleSeq.rna2dna(miR.getSeq()));
                            }
                        }
                    brCounts.close();
 
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
        
        // write uniqMirNames
        logger.info("writing uniq isomiR list");
        try{
            BufferedWriter bwUI = new BufferedWriter(new FileWriter(new File(uniqIsomiRsFilePath.toString())));
              for(String key: uniqMirNames.keySet()){
                bwUI.write(key + "\t"  + uniqMirNames.get(key)+ System.getProperty("line.separator") );
              }
            bwUI.close();
        }
        catch(IOException exIO){
          exIO.printStackTrace();
          logger.info("error writing unique isomiR list File <" + uniqIsomiRsFilePath.toString() + ">\n" + exIO);
          throw new IOException(STEP_ID_STRING + "error writing isomiR list File <" + uniqIsomiRsFilePath.toString() + ">");
        }
        
        logger.info(STEP_ID_STRING + ": completed");
 
    }
    
    
    
    /**
     * Verify Input Data for parsing SAM file for miRNAs
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
   * @return the groupBySeedRegion
   */
  public Boolean getGroupBySeedRegion() {
    return groupBySeedRegion;
  }

  /**
   * @param groupBySeedRegion the groupBySeedRegion to set
   */
  public void setGroupBySeedRegion(Boolean groupBySeedRegion) {
    this.groupBySeedRegion = groupBySeedRegion;
  }


}
