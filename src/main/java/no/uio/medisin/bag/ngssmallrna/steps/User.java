/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package no.uio.medisin.bag.ngssmallrna.steps;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author sr
 */
public class User {
 
/*
    
    private int age = 29;
	private String name = "mkyong";
	private List<String> messages = new ArrayList<String>() {
		{
			add("msg 1");
			add("msg 2");
			add("msg 3");
		}
	};
    */
    private int age;
    private String name;
    private List<String> messages = new ArrayList<String>();

    private String name2;
    private List <NGSRunStepJSONData> steps = new ArrayList<>();
        
 
	//getter and setter methods
 
	@Override
	public String toString() {
		return "User [age=" + age + ", name=" + name + ", " +
				"messages=" + messages + "]";
	}
}
