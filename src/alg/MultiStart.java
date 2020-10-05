package alg;


import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.ThreadLocalRandom;

import umontreal.iro.lecuyer.probdist.Distribution;
import umontreal.iro.lecuyer.probdist.LognormalDist;
import umontreal.iro.lecuyer.randvar.LognormalGen;
import umontreal.iro.lecuyer.randvar.RandomVariateGen;
import umontreal.iro.lecuyer.rng.LFSR113;
import umontreal.iro.lecuyer.rng.RandomStream;
import umontreal.iro.lecuyer.rng.RandomStreamBase;

/**
 * Iteratively calls the RandCWS and saves the best solution.
 * @author Angel A. Juan - ajuanp(@)gmail.com
 * @version 13810
 */
public class MultiStart 
{
    /* 0. Instance fields and class constructor */
    private Test aTest;
    private Inputs inputs; // contains savingsList too
    private Random rng;
    private Solution[] initialSols = new Solution[6]; // (0):0% (1):25% ... (5):decent.
    private Solution initialSol = null;
    private Solution bestSol = null;
    private Solution bestSolDist = null;
    private Solution bestSolVio = null;
    private Solution newSol = null;
    private Solution baseSol = null;
    private Outputs outputs = new Outputs();
    RandomStreamBase stream = new LFSR113(); // L'Ecuyer stream
    Node[] auxTempNodes; 
    List<Double> inventoryCost= new ArrayList<>();
    List<Double> routingCost= new ArrayList<>();
    ArrayList<Solution> outList = new ArrayList<Solution>();
    
    MultiStart(Test myTest, Inputs myInputs, Random myRng)
    {
        aTest = myTest;
        inputs = myInputs;
        rng = myRng;
    }
    
    
    

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%% Deterministic %%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    public Solution multiCW(Solution sol, double alpha){
    	Solution newsol = null;
        long start = ElapsedTime.systemTime();
        double elapsed = 0.0;
        boolean firstime = true;
        
  	    initialSol = new Solution(sol); //Lo llamo 1 vez fuera del bucle para generar una solución inicial
  	    System.out.println("Solucion inicial: " + initialSol.getTotalScore());
  	    baseSol = initialSol;
        bestSol = initialSol;  
        InputsManager.generateSavingsList(inputs,alpha);

   
        while( elapsed < aTest.getMaxTime() )
        {     
    	 newSol = new Solution(CWS.solve(inputs,aTest,rng,1));  //APlico C&W bias
    	 Collections.sort(newSol.getRoutes());  //Me quedo con las N mejores rutas (igual al numero de vehiculos)
		 double totalCost = 0.0;
		 double totalProfit = 0.0;
		 int used_veh = Math.min(inputs.getVehNumber(),newSol.getRoutes().size());
		 newSol.sliceSolution(used_veh);
		 for(int i = 0; i < used_veh; i++){
				Route r = newSol.getRoutes().get(i);
				totalCost += r.getCosts();
				totalProfit += r.getScore();
		}	
		 newSol.setTotalCosts(totalCost);
		 newSol.setTotalScore(totalProfit);
 
    	 if(newSol.getTotalScore() > baseSol.getTotalScore()){
    		 baseSol = newSol;
    		 
    		 if(newSol.getTotalScore() > bestSol.getTotalScore()){
    			 bestSol = newSol;
    			 bestSol.setTime(elapsed);
    			 System.out.println("Mejoro SOL : " + bestSol.getTotalScore());
    		 }
    		 
    	 }

          elapsed = ElapsedTime.calcElapsed(start, ElapsedTime.systemTime());
        }
    	    	 
    	 return bestSol;
    }
    
    
    
   

    
    private static boolean merge(Route route, Edge edge, Node node, Test aTest, Inputs inputs, double costLastNode) {
		
    	boolean ismerge = false;
    	
		Node p = edge.getOrigin();
		Node q = edge.getEnd();
		String key1=p.getId()+","+ node.getId();
		String key2=node.getId()+","+q.getId();
		String key3=p.getId()+","+q.getId();
		Edge e1=inputs.getEdgesDirectory().get(key1);
		
		//if(e1==null) {
			System.out.println("e1 " +e1.toString());
		//}
		Edge e2=inputs.getEdgesDirectory().get(key2);
		//if(e2==null) {
			System.out.println("e2 "+e2.toString());
		//}
		Edge e3=inputs.getEdgesDirectory().get(key3);
		//if(e3==null) {
			System.out.println("e3 "+e3.toString());
		//}
		
		double cost = (e1.getCosts() +
				e2.getCosts() - 
				e3.getCosts() + 0);
		if( (route.getCosts() + cost + costLastNode ) <= inputs.gettMax()){
			ismerge = true;
			route.setCosts((route.getCosts() + cost));
		}

		return ismerge;
	}
  
    
    
    public static int getRandomPosition(double alpha, double beta, Random r,int size) {
		 double randomValue = alpha + (beta - alpha) * r.nextDouble();
		 int index = (int) (Math.log(r.nextDouble()) / Math.log(1 - randomValue));
		 index = index % size;
		 return index;
	}
    
    
    
    
	public static int getRandomPosition( double beta, Random r, int size) {
		 int index = (int) (Math.log(r.nextDouble()) / Math.log(1 - beta));
         index = index % size;
     return index;
	}
    

	
	
	
public static TreeSet<PairBest> calEvaluationFunction(Set<Node> nodes,Edge edge, Route r, Inputs inputs) {
		
		TreeSet<PairBest> distanceNodes = new TreeSet<PairBest>();
		
		Node p = edge.getOrigin();
		Node q = edge.getEnd();
		for(Node inode: nodes){
			String key1=p.getId()+","+inode.getId();
			Edge e1=inputs.getEdgesDirectory().get(key1);
			String key2=inode.getId()+","+q.getId();
			Edge e2=inputs.getEdgesDirectory().get(key2);
			String key3=p.getId()+","+q.getId();
			Edge e3=inputs.getEdgesDirectory().get(key3);
			
			double cost = (e1.getCosts() + e2.getCosts() - e3.getCosts() + 0) / inode.getProfit();
			PairBest nodeDist = new PairBest(cost,inode);
			distanceNodes.add(nodeDist);
		}

		return distanceNodes;
	}
	
	

public static TreeSet<PairBestProfit> calEvaluationFunctionProfit(Set<Node> nodes,Edge edge, Route r, Inputs inputs) {
	
	TreeSet<PairBestProfit> distanceNodes = new TreeSet<PairBestProfit>();
	
	Node p = edge.getOrigin();
	Node q = edge.getEnd();
	for(Node inode: nodes){
		//double cost = (edge.calcCosts(p, inode) + edge.calcCosts(inode, q) - edge.calcCosts(p, q) + 0) / inode.getProfit();
		double cost = inode.getProfit();
		PairBestProfit nodeDist = new PairBestProfit(cost,inode);
		distanceNodes.add(nodeDist);
	}

	return distanceNodes;
}


   
    
	public static PairBest obtainDistanceNode(int position, TreeSet<PairBest> tree) {
		
		Iterator it = tree.iterator();
        int k = 0;
        PairBest p = null;
        
		while(it.hasNext()) {
            p = (PairBest) it.next();
            if(position == k){
            	return p;
            }
        k++;
		}    
		return p;
	}
	
	
	
	public static PairBestProfit obtainProfitNode(int position, TreeSet<PairBestProfit> tree) {
		
		Iterator it = tree.iterator();
        int k = 0;
        PairBestProfit p = null;
        
		while(it.hasNext()) {
            p = (PairBestProfit) it.next();
            if(position == k){
            	return p;
            }
        k++;
		}    
		return p;
	}
	
	
	
	
	
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%% Heuristica C&W local %%%%%%%%%%%%%%%%%%%%%%
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	public Solution solveCWS(){
		InputsManager.generateDepotEdges(inputs);
		Solution best = new Solution();
		best.setTotalCosts(Double.MAX_VALUE);
		
		
		///1- Aplicaciones heurísticaC&W
		double bAlpha = 0.0;
		for (double alpha = 0; alpha <= 1; alpha+=0.1){
			InputsManager.generateSavingsList(inputs,alpha);
			Solution detSol = CWS.solve(inputs,aTest,rng,0);
			Collections.sort(detSol.getRoutes());
			double totalCost = 0.0;
			double totalProfit = 0.0;
			int used_veh = Math.min(inputs.getVehNumber(),detSol.getRoutes().size());

			detSol.sliceSolution(used_veh);

			for(int i = 0; i < used_veh; i++){
				Route r = detSol.getRoutes().get(i);
				totalCost += r.getCosts();
				totalProfit += r.getScore();
			}
			
			detSol.setTotalCosts(totalCost);
			detSol.setTotalScore(totalProfit);
			
			if(totalProfit > best.getTotalScore() || 
			   (totalProfit == best.getTotalScore() && totalCost < best.getTotalCosts())){
				best = detSol;
				bAlpha = alpha;
			}
			
		}
		
		
		
		///2- Aplicaciones heurísticaC&W proceso constructivo para la mejor alpha
		best = multiCW(best,bAlpha); //Ejecución determinista			
		return best;
	}
	

    
    
  
}










