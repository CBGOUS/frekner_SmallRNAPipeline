/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.pipeline;

import junit.framework.Assert;
import org.testng.annotations.AfterClass;
import org.testng.annotations.AfterMethod;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

/**
 * 
 * @author simonray
 */
public class SampleDataEntryNGTest {
    
    public SampleDataEntryNGTest() {
    }

    
    @Test(expectedExceptions = Exception.class)
    public void missingFastq1() {
        SampleDataEntry sed = new SampleDataEntry(null, "fastq2", "source", "condition", "time", "note");
            
    }

    @Test
    public void missingFastq2() {
        SampleDataEntry sed = new SampleDataEntry("fastq1", null, "source", "condition", "time", "note");
        //Assert.assertNotNull(sed.getFastqFile2());
    }
    
        
    @Test
    public void missingSource()
    {
        SampleDataEntry sed = new SampleDataEntry("fastq1", null, "source", "condition", "time", "note");     
        Assert.assertNotNull(sed.getDataSource());
    }   
    
    @Test
    public void missingCondition()
    {
        SampleDataEntry sed = new SampleDataEntry("fastq1", null, "source", "condition", "time", "note");     
        Assert.assertNotNull(sed.getCondition());
    }   
    
    @Test
    public void missingTime()
    {
        SampleDataEntry sed = new SampleDataEntry("fastq1", null, "source", "condition", "time", "note");     
        Assert.assertNotNull(sed.getTime());
    }   
    
    @Test
    public void missingNote()
    {
        SampleDataEntry sed = new SampleDataEntry("fastq1", null, "source", "condition", "time", "note");     
        Assert.assertNotNull(sed.getNote());
    }   
    
    
    @BeforeClass
    public static void setUpClass() throws Exception {
    }

    @AfterClass
    public static void tearDownClass() throws Exception {
    }

    @BeforeMethod
    public void setUpMethod() throws Exception {
    }

    @AfterMethod
    public void tearDownMethod() throws Exception {
    }

    /**
     * Test of toString method, of class SampleDataEntry.
     */
    @Test
    public void testToString() {
    }

    /**
     * Test of getFastqFile1 method, of class SampleDataEntry.
     */
    @Test
    public void testGetFastqFile1() {
    }

    /**
     * Test of setDataFile method, of class SampleDataEntry.
     */
    @Test
    public void testSetDataFile() {
    }

    /**
     * Test of getDataSource method, of class SampleDataEntry.
     */
    @Test
    public void testGetDataSource() {
    }

    /**
     * Test of setDataSource method, of class SampleDataEntry.
     */
    @Test
    public void testSetDataSource() {
    }

    /**
     * Test of getCondition method, of class SampleDataEntry.
     */
    @Test
    public void testGetCondition() {
    }

    /**
     * Test of setCondition method, of class SampleDataEntry.
     */
    @Test
    public void testSetCondition() {
    }

    /**
     * Test of getTime method, of class SampleDataEntry.
     */
    @Test
    public void testGetTime() {
    }

    /**
     * Test of setTime method, of class SampleDataEntry.
     */
    @Test
    public void testSetTime() {
    }

    /**
     * Test of getNote method, of class SampleDataEntry.
     */
    @Test
    public void testGetNote() {
    }

    /**
     * Test of setNote method, of class SampleDataEntry.
     */
    @Test
    public void testSetNote() {
    }

    /**
     * Test of getFastqFile2 method, of class SampleDataEntry.
     */
    @Test
    public void testGetFastqFile2() {
    }

    /**
     * Test of setFastqFile2 method, of class SampleDataEntry.
     */
    @Test
    public void testSetFastqFile2() {
    }
    
}
