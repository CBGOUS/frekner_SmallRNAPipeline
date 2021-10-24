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
import no.uio.medisin.bag.core.mapping.SAMEntry;
import no.uio.medisin.bag.core.mathstuff.FWHM;

import no.uio.medisin.bag.core.mirna.MiRNA;
import no.uio.medisin.bag.core.mirna.PreMiRNA;
import no.uio.medisin.bag.core.mirna.PreMiRNASet;
import no.uio.medisin.bag.core.sequence.Strand;
import no.uio.medisin.bag.ngssmallrna.pipeline.ReferenceDataLocations;
import no.uio.medisin.bag.ngssmallrna.pipeline.SampleDataEntry;
import static no.uio.medisin.bag.ngssmallrna.steps.NGSStep.FILESEPARATOR;
import org.apache.commons.lang3.StringUtils;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;


/**
 *   1. parse SAM file to extract and process the miRNA reads associated with a 
 *      set of pre-miRNAs
 * 
 *   Input is a SAM file
 * 
 * @author sr
 */

public class StepParseSAMForPreMiRNAReads extends NGSStep{
    
    static Logger                   logger                          = LogManager.getLogger();
    
    public static final NGSStepSubclass STEP_SUBCLASS               = NGSStepSubclass.DATABLE;
    
    public  static final String     STEP_ID_STRING                  = "ParseSAMForPreMiRNAReads";
    private static final String     ID_BLEED                        = "bleed";
    private static final String     ID_ISOMIRS                      = "analyzeIsomirs";
    private static final String     ID_MIRBASE_VERSION              = "mirbaseVersion";
    private static final String     ID_REF_GENOME                   = "host";
    private static final String     ID_BASELINE                     = "baselinePercent";
        

    private static final String     RAW_INPUT_EXTENSION             = ".fastq.gz";
    private static final String     INFILE_EXTENSION                = ".trim.clp.gen.sam";
    private static final String     PREMIR_DETAILS_EXTENSION        = ".trim.clp.gen.premir_details.txt";
    private static final String     PREMIRNASETS_EXTENSION          = ".trim.clp.gen.mircounts.tsv";
    
    public  static final String     PREMIRNA_START_ENTRY            = "<START_PREMIRNA>";
    public  static final String     PREMIRNA_END_ENTRY              = "<END_PREMIRNA>";
    public  static final String     PREMIRNA_NAME_SECTION           = "PREMIRNA_NAME:";
    public  static final String     PREMIRNA_FWHM_SECTION           = "FWHM:";
    public  static final String     PREMIRNA_LENGTH_SECTION         = "LENGTH_DISTRIBUTION:";
    public  static final String     PREMIRNA_START_SECTION          = "START_DISTRIBUTION:";
    public  static final String     PREMIRNA_LEN_FWHM               = "LEN_FWHM:";
    public  static final String     PREMIRNA_START_FWHM             = "START_FWHM:";
    public  static final String     PREMIRNA_UPPERLEN_FWHM          = "UPPER_LEN_FWHM:";
    public  static final String     PREMIRNA_UPPERSTART_FWHM        = "UPPER_START_FWHM:";
    public  static final String     PREMIRNA_LOWERLEN_FWHM          = "LOWER_LEN_FWHM:";
    public  static final String     PREMIRNA_LOWERSTART_FWHM        = "LOWER_START_FWHM:";
    


    private PreMiRNASet             preMiRNAHitList;
    PreMiRNASet                     mirBasePreMiRNASet              = new PreMiRNASet();           
    
    private int                     locationBleed                   = 2;
    private Boolean                 generateMetrics                  = false;
    private int                     miRBaseRelease                  = 20;
    private String                  referenceGenome                 = "";
    private int                     baselinePercent                 = 5;
    



    public StepParseSAMForPreMiRNAReads(){
        classSubtype = NGSStepSubclass.DATABLE;
    }
    
    /**
     * 
     * @param sid StepInputData
     * 
     */
    public StepParseSAMForPreMiRNAReads(InputDataForStep sid){
        classSubtype = NGSStepSubclass.DATABLE;
        stepInputData = sid;
    }
    
    @Override
    public void parseStepParameters() throws Exception{
      
    }
    
    
    
    @Override
    public String shortStepDescription(){
      return "parse SAM file to extract and process the miRNA reads associated with a set of pre-miRNAs";
    }
    
    
    
    @Override
    public String longStepDescription(){
      return "parse SAM file to extract and process the miRNA reads associated with a set of pre-miRNAs\n"
              + "Input is a SAM file.\n";
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
        chk = checkParameter("Integer", ID_MIRBASE_VERSION, Integer.toString((Integer)configData.get(ID_MIRBASE_VERSION)), "0", "NA", logger);
        if(chk!=null)
            this.setMiRBaseRelease((Integer) configData.get(ID_MIRBASE_VERSION));

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
    
        String miRBaseGFFFile = this.cleanPath(getStepInputData().getDataLocations().getMirbaseFolder() 
                + FILESEPARATOR + this.getMiRBaseRelease() + FILESEPARATOR + this.getReferenceGenome() + ".gff3");
        
        Boolean fA = new File(outFolder).mkdir();       
        if (fA) logger.info("created output folder <" + outFolder + "> for results" );
        String samLine = null;
        String samInputFile = "";

        PreMiRNASet mirGeneDBList = new PreMiRNASet();
        String mirgdbLine = "";
        String mirGeneDBFile = this.cleanPath(getStepInputData().getDataLocations().getMirgenedbfolder()
            + FILESEPARATOR + this.getReferenceGenome() + ".mirgenedb.tsv");
        try{
        BufferedReader brMG = new BufferedReader(new FileReader(new File(mirGeneDBFile)));
            brMG.readLine(); // header line
            while((mirgdbLine=brMG.readLine())!=null){
                PreMiRNA preMiRNA = new PreMiRNA();
                preMiRNA.parseMirGeneDBEntry(mirgdbLine);
                mirGeneDBList.addPriMiRNA(preMiRNA);
            }
        brMG.close();
        }
        catch(IOException exIO){
            logger.info("Exception parsing MirGeneDB file <" + mirGeneDBFile  + ">");
            logger.info(exIO);
            throw new IOException("Exception parsing MirGeneDB file <" + mirGeneDBFile  + ">");
        }

        String genomeFastaFile = this.cleanPath(
                getStepInputData().getDataLocations().getGenomeRootFolder()
                        + FILESEPARATOR + this.getReferenceGenome() 
                        + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_PATH 
                        + FILESEPARATOR + "genome.fa");

        mirBasePreMiRNASet.loadMiRBaseData(this.getReferenceGenome(), miRBaseGFFFile, genomeFastaFile);
        
        /*
            match the miRBase entry to the miRGeneDB entry
        */
        ArrayList<String> preMiRNAResults = new ArrayList<>();
        preMiRNAResults.add("miRBase Name|miRBase ID|miRGeneDB Name|miRGeneDB Gene Node\t");
        for(PreMiRNA preMiRNA: mirBasePreMiRNASet.getPreMiRNAList()){
            Boolean miRNAmatch = false;
            for(MiRNA miRNA: preMiRNA.getMiRNAList()){
                PreMiRNA mirGeneDBHit = mirGeneDBList.doesListContainMiRNA(miRNA.getMimatID());
                if(mirGeneDBHit!=null){
                    preMiRNAResults.add(preMiRNA.getName() + "|" + preMiRNA.getMiID()
                            + "|" + mirGeneDBHit.getMirGeneDBGeneName() 
                            + "|" + mirGeneDBHit.getMirGeneDBNodeOfOriginGene() + "\t");
                    miRNAmatch = true;
                    break;
                } 
                
            }
            if(!miRNAmatch){
                preMiRNAResults.add(preMiRNA.getName() + "|" + preMiRNA.getMiID() 
                            + "|||\t");
            }
        }
        
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
                long totalCounts = 0;
                samLine = null;
                BufferedReader brSAM = new BufferedReader(new FileReader(new File(samInputFile)));

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
                    SAMEntry samEntry = new SAMEntry(samLine);

                    if(samEntry.isHeaderLine()) continue;

                    if(samLine.split("\t")[SAMEntry.COL_QNAME].contains("-")){
                        totalCounts += Integer.parseInt(samLine.split("\t")[SAMEntry.COL_QNAME].split("-")[1]);                            
                    }else{
                        totalCounts++;
                    }

                    if(samEntry.isMappedRead()){

                        if (samEntry.getStrand()==Strand.MINUS) 
                            preMatchCount3++;
                        else
                            preMatchCount5++;

                        PreMiRNA preMiRBaseHit = mirBasePreMiRNASet.doesReadOverlapKnownPreMiRNA(samEntry.getStartPos(), samEntry.getEndPos(), 
                                samEntry.getrName(), samEntry.getStrand().toString(), bleed);

                        if (preMiRBaseHit != null){

                            if (samEntry.getStrand()==Strand.MINUS) 
                                matchCount3 += Integer.parseInt(samLine.split("\t")[SAMEntry.COL_QNAME].split("-")[1]);
                            else
                                matchCount5 += Integer.parseInt(samLine.split("\t")[SAMEntry.COL_QNAME].split("-")[1]);
                            if(samEntry.getSeq().equals("")){
                                logger.info("eh?");
                            }
                            if(Integer.parseInt(samEntry.getqName().split("-")[1]) > bleed){
                                MiRNA miRRead = new MiRNA(samEntry.getqName(), 
                                        samEntry.getrName(), samEntry.getStartPos()+1, samEntry.getEndPos()+1, 
                                        samEntry.getStrand().toString(), preMiRBaseHit.getMiID() + ":"  + samEntry.getqName(), 
                                        preMiRBaseHit.getMiID(), samEntry.getSeq());
                                miRRead.setReadCount(Integer.parseInt(samEntry.getqName().split("-")[1]));

                                mirBasePreMiRNASet.addMiRNAtoPreMiRNA(miRRead, preMiRBaseHit.getMiID());
                                logger.debug("add: " + miRRead.getName() + " to " + preMiRBaseHit.getMiID() + "[" + preMiRBaseHit.getStrand() + "]( " 
                                        + mirBasePreMiRNASet.getEntryByMiID(preMiRBaseHit.getMiID()).getReadCount() + ")");
                                logger.debug("have " + mirBasePreMiRNASet.getPreMiRNAList().size() + " entries");
                            }

                        }

                    }
                }
                logger.info("  total mapped counts = " + totalCounts);
                
                /*
                    the following is rather approximate.
                    apparently, for 5,000,000 reads, the lowest detectable by qPCR is 50. so, we divide total counts by 100000.
                    there has to be a better way....
                */
                Double minCountsForSingleFeature = (double) totalCounts /100000.0; 
                logger.info((matchCount5 + matchCount3) + " reads (" + matchCount5 + " 5'" + "/" + matchCount3 + " 3' ) were mapped to pre-miRNAs");

                if(generateMetrics){
                    logger.info("  calculate SAM metrics");
                    totalCounts = mirBasePreMiRNASet.countUpReadsInAllPreMiRNAs();
                    logger.info(totalCounts);


                    int p=0;
                    String countString = sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " Read Counts" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " miR Counts" + "\t";
                    String plusStartString = sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 5' min read start" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 5' max read start" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 5' start average" + "\t";
                    String minusStartString = sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 3' min read start" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + "3' max read start" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + "3' start average" + "\t";
                    String plusLenString = sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 5' min read length" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 5' max read length" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 5' average read length" + "\t";
                    String minusLenString = sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 3' min read length" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 3' max read length" + "\t"
                            + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, "") + " 3' average read length" + "\t";
                    preMiRNAResults.set(p, preMiRNAResults.get(p).concat(countString 
                            + plusStartString 
                            + minusStartString
                            + plusLenString
                            + minusLenString
                            ));
                    

                    p++;
                    
                    
                    logger.info("  write pre-miRNA details");
                    try{
                        String  preMiRNADetailsFile = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, PREMIR_DETAILS_EXTENSION);
                        logger.info("write to " + preMiRNADetailsFile);
                        BufferedWriter brPreDetails = new BufferedWriter(new FileWriter(new File(preMiRNADetailsFile)));
                            for(PreMiRNA preMiRHit: mirBasePreMiRNASet.getPreMiRNAList()){
                                brPreDetails.write("\n\n+" + StringUtils.repeat("-", 80) + "+\n\n");
                                brPreDetails.write(PREMIRNA_START_ENTRY + "\n");
                                brPreDetails.write(PREMIRNA_NAME_SECTION + "\t" + preMiRHit.getMiID() + "|" + preMiRHit.getName() + "\n");
                                logger.debug(preMiRHit.getName());
                                logger.debug("\n" + preMiRHit.prettyPlotReadsWithSeq() + "\n");
                                double lenFWHM = FWHM.firstMinusLastReadDistance(preMiRHit.detailLengthRange(1, 40), 40);
                                logger.debug(lenFWHM);
                                double startFWHM = FWHM.firstMinusLastReadDistance(preMiRHit.detailStartRange(0, 99), preMiRHit.getUpperEnd()-preMiRHit.getUpperStart()+1);
                                logger.debug(startFWHM);
                                double upperLenFWHM = FWHM.firstMinusLastReadDistance(preMiRHit.detailLengthRangeFiveArm(1, 40), 40);
                                logger.debug(upperLenFWHM);
                                double upperStartFWHM = FWHM.firstMinusLastReadDistance(preMiRHit.detailStartRangeFiveArm(0, 99), preMiRHit.getUpperEnd()-preMiRHit.getUpperStart()+1);
                                logger.debug(upperStartFWHM);
                                double lowerLenFWHM = FWHM.firstMinusLastReadDistance(preMiRHit.detailLengthRangeThreeArm(1, 40), 40);
                                logger.debug(lowerLenFWHM);
                                double lowerStartFWHM = FWHM.firstMinusLastReadDistance(preMiRHit.detailStartRangeThreeArm(0, 99), preMiRHit.getUpperEnd()-preMiRHit.getUpperStart()+1);
                                logger.debug(lowerStartFWHM);
                                brPreDetails.write(PREMIRNA_FWHM_SECTION + "\n");
                                brPreDetails.write(PREMIRNA_LEN_FWHM        + "\t" + lenFWHM + "\n");
                                brPreDetails.write(PREMIRNA_START_FWHM      + "\t" + startFWHM + "\n");
                                brPreDetails.write(PREMIRNA_UPPERLEN_FWHM   + "\t" + upperLenFWHM + "\n");
                                brPreDetails.write(PREMIRNA_UPPERSTART_FWHM + "\t" + upperStartFWHM + "\n");
                                brPreDetails.write(PREMIRNA_LOWERLEN_FWHM   + "\t" + lowerLenFWHM + "\n");
                                brPreDetails.write(PREMIRNA_LOWERSTART_FWHM + "\t" + lowerStartFWHM + "\n");


                                brPreDetails.write("\n" + preMiRHit.prettyPlotReadsWithSeq() + "\n");

                                brPreDetails.write("\n" + PREMIRNA_LENGTH_SECTION + "\n");
                                brPreDetails.write(preMiRHit.detailLengthRange(1, 40).toString());

                                brPreDetails.write("\n" + PREMIRNA_START_SECTION + "\n");
                                brPreDetails.write(preMiRHit.detailStartRange(0, 99).toString());
                                brPreDetails.write("\n\n+" + StringUtils.repeat("-", 80) + "+\n\n");
                                brPreDetails.write(PREMIRNA_END_ENTRY + "\n");

                                preMiRNAResults.set(p, preMiRNAResults.get(p) + preMiRHit.getReadCount() + "\t" + preMiRHit.getMiRNAList().size() + "\t"
                                        + startFWHM + "\t"
                                        + lenFWHM  + "\t");
                                p++;
                            }
                        brPreDetails.close();
                    }
                    catch(IOException ioEX){
                        
                    }
                }

                
                logger.info("  write pre-miRNA/miRNA sets");

                String  preMiRNASetsFile  = outFolder + FILESEPARATOR + sampleData.getFastqFile1().replace(RAW_INPUT_EXTENSION, PREMIRNASETS_EXTENSION);

                BufferedWriter brPriSets  = new BufferedWriter(new FileWriter(new File(preMiRNASetsFile)));
                    for(PreMiRNA preMiRHit: mirBasePreMiRNASet.getPreMiRNAList()){
                        brPriSets.write(preMiRHit.getMiID() + "\n");
                        for( MiRNA miRHit: preMiRHit.getMiRNAList()){
                            brPriSets.write("\t" + miRHit.getMimatID() + "\t" + miRHit.getName() + "\t"  
                                    + miRHit.getMiStart() + "\t" + miRHit.getMiEnd() 
                                    + "\t" + miRHit.getStrand() + "\t" + miRHit.getReadCount() + "\n");
                        }
                    }
                brPriSets.close();
                    
                mirBasePreMiRNASet.removeNGSReadEntries();

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
        
        if(generateMetrics){
            
            String allCountsFile   = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".readsandmir.counts.tsv";
            //String summaryFile      = outFolder + FILESEPARATOR + getStepInputData().getProjectID() + ".disp.summary.tsv";
            logger.info("write dispersions to file <" + allCountsFile + ">");
            try{
                BufferedWriter bwAc = new BufferedWriter(new FileWriter(new File(allCountsFile)));
                //BufferedWriter bwSm = new BufferedWriter(new FileWriter(new File(summaryFile)));  
                /*
                    here we want to add in the miRGene information
                    we can only do this by matching the MIMATID, we have nothing else to go on
                    The problem is made more complicated because a miRNA can be degenerate, so we have
                    to use additional information to figure out which miRGeneDB entry is correct
                */
                for(String preMiRNAresultLine: preMiRNAResults){
                    
                    bwAc.write(preMiRNAresultLine + "\n");
                }
                //bwSm.close();
                bwAc.close();
            }
            catch(IOException exIO){
                logger.info("error writing isomiR dispersion File <" + allCountsFile + ">\n" + exIO);
                throw new IOException(STEP_ID_STRING + "error writing isomiR dispersion File <" + allCountsFile + ">");
            }
           
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
                
        String genomeFastaFile = this.cleanPath(
                getStepInputData().getDataLocations().getGenomeRootFolder()
                        + FILESEPARATOR + this.getReferenceGenome() 
                        + FILESEPARATOR + ReferenceDataLocations.ID_REL_WHOLE_GENSEQ_PATH 
                        + FILESEPARATOR + "genome.fa");
        if (new File(genomeFastaFile).exists()==false){
            logger.error("no genome fasta reference file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + genomeFastaFile + ">");
            throw new IOException("no genome fasta reference file was found for mirBase HOST:<" 
                    + this.getReferenceGenome() + "> VERSION: <"+ this.getMiRBaseRelease() + "> at location <" 
                    + genomeFastaFile + ">");
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
        return generateMetrics;
    }

    /**
     * @param analyzeIsomirs the analyzeIsomirs to set
     */
    public void setAnalyzeIsomirs(Boolean analyzeIsomirs) {
        this.generateMetrics = analyzeIsomirs;
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

}
