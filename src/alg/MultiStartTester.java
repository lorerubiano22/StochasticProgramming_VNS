package alg;
import java.awt.List;
import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;

import umontreal.iro.lecuyer.rng.LFSR113;
import umontreal.iro.lecuyer.rng.RandomStreamBase;

import java.io.IOException;
import java.io.PrintWriter;

/**
 * @author Angel A. Juan - ajuanp(@)gmail.com
 * @version 130112
 */
public class MultiStartTester
{
	final static String inputFolder = "inputs";
	final static String outputFolder = "outputs";
	final static String testFolder = "tests";
	final static String fileNameTest = "test2Run.txt"; 
	final static String sufixFileNodes = ".txt";
	final static String sufixFileVehicules = "_input_vehicles.txt";
	final static String sufixFileOutput = "_outputs.txt";



	public static void main( String[] args )
	{
		BKS.readBKS("BKS.csv");
		ArrayList<Outputs> outList = new ArrayList<Outputs>();
		ArrayList<Outputs> outListDist = new ArrayList<Outputs>();

		ArrayList<ArrayList<Outputs>> solToPrint = new ArrayList<ArrayList<Outputs>>(2);


		System.out.println("****  WELCOME TO THIS PROGRAM  ****");
		long programStart = ElapsedTime.systemTime();

		/* 1. GET THE LIST OF TESTS TO RUN FORM "test2run.txt"
              aTest = instanceName + testParameters */
		String testsFilePath = testFolder + File.separator + fileNameTest;
		ArrayList<Test> testsList = TestsManager.getTestsList(testsFilePath);

		/* 2. FOR EACH TEST (instanceName + testParameters) IN THE LIST... */
		int nTests = testsList.size();
		System.out.println("Number of tests to run: " + nTests);
		Test aTest = null;
		int num_opt = 0;
		int num_sol = 0;
		for( int k = 0; k < nTests; k++ )
		{   
			aTest = testsList.get(k);
			//System.out.println("\n# STARTING TEST " + (k + 1) + " OF " + nTests);
			System.out.print("Start test for Instance: " + aTest.getInstanceName());
			System.out.println();

			// 2.1 GET THE INSTANCE INPUTS (DATA ON NODES AND VEHICLES)
			// "instanceName.txt" contains data on nodes
			String inputNodesPath = inputFolder + File.separator + aTest.getInstanceName() + sufixFileNodes;


			// Read inputs files (nodes) and construct the inputs object
			Inputs inputs = InputsManager.readInputs(inputNodesPath);

			//Tengo las distancias del depot de salida a los customers odenadas de mayor a menor
			//InputsManager.generateDepotEdges(inputs);

			//USE THE MULTI-START ALGORITHM TO SOLVE THE INSTANCE
			Random rng = new Random(aTest.getSeed());
			RandomStreamBase stream = new LFSR113(); // L'Ecuyer stream
			aTest.setRandomStream(stream);

			VNS alg2 = new VNS(aTest,inputs,rng);
			Solution bestSolution = null;
			Outputs output = null;



			bestSolution = alg2.solveDet();
			output = new Outputs();

			output.setOBSol(bestSolution);
			output.setInstanceName(aTest.getInstanceName());
			output.setseed(aTest.getSeed());
			outList.add(output);
			Integer best_score = BKS.bestSolution(aTest.getInstanceName());  
		}

		Outputs.printSolST(outList,aTest);

		Outputs.printSol(outList);



		/* 3. END OF PROGRAM */
		System.out.println("****  END OF PROGRAM, CHECK OUTPUTS FILES  ****");
		long programEnd = ElapsedTime.systemTime();
		System.out.println("Total elapsed time = "
				+ ElapsedTime.calcElapsedHMS(programStart, programEnd));
	}   

}