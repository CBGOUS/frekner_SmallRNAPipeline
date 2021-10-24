/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import java.util.ArrayList;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

/**
 * stores a set of Conditions, and implements specific operations on the set
 * 
 * @author sr
 */
public class ConditionSet {
    
    Logger logger = LogManager.getRootLogger();
    private ArrayList<ConditionEntry>     pairedList = new ArrayList();

    
    
    public void addEntry(String n, String d1){
        pairedList.add(new ConditionEntry(n, d1, ""));
    }
    
    
    public void updateEntry(String n, String d2){
        for(ConditionEntry conditionEntry: getPairedList()){
            if (conditionEntry.getName().equals(n)){
                conditionEntry.setDataset2(d2);
            }
        }
    }

    /**
     * remove any entries that don't have both datasets defined
     * 
     * @return 
     */
    public int removeUnpairedEntries(){
        
        int startSize = pairedList.size();
        for(int i=0; i<pairedList.size(); i++){
            ConditionEntry conditionEntry = (ConditionEntry) pairedList.get(i);
            if (conditionEntry.getDataset1().isEmpty() || conditionEntry.getDataset2().isEmpty()){
                logger.info("remove unpaired entry " + conditionEntry.getName());
                pairedList.remove(i--);
            }
        }
        return startSize - pairedList.size();
    }
    /**
     * @return the pairedList
     */
    public ArrayList<ConditionEntry> getPairedList() {
        return pairedList;
    }
    
    /**
     * return how many paired entries in the set
     * 
     * @return 
     */
    public int getNumberOfEntries(){
        return pairedList.size();
    }
    
    public String printPairedList(){
        String returnString = "";
        for(ConditionEntry conditionEntry: pairedList){
            returnString = returnString.concat(conditionEntry.printEntry() + "\n");
        }
        return returnString;
    }
    
    
}
