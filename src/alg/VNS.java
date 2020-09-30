package alg;

import java.util.*;

import umontreal.iro.lecuyer.randvar.LognormalGen;
import umontreal.iro.lecuyer.rng.RandomStream;

/**
 * Created by lluc on 2/06/17.
 * This class encapsulates the metaheuristic VNS
 */
public class VNS {

	private static final int P_MIN = 1 ;
	private static final int P_MAX = 100 ;
	private static final int P_STEP = 2 ;
	private static Test aTest;
	private static Inputs inputs;
	private Random rng;
	private double bestAlpha;
	private Solution initialSolution;
	private Solution bestSolution;
	private static final double EPSILON = 10E-6;
	RouteCache  routecache = new RouteCache();
	private double[][] cumulativeObjectivefunction=null;



	VNS(Test myTest, Inputs myInputs, Random myrng)
	{
		aTest = myTest;
		inputs = myInputs;
		rng = myrng;
		cumulativeObjectivefunction= new double[aTest.getLongSim()][2];
	}





	/**
	 * It applies VNS metaheuristic to solve the problem specified at attribute inputs, using test instance class aTest
	 * @return The best solution found by VNS
	 */
	public Solution solveDet(){
		for(int scenario=0;scenario<aTest.getLongSim();scenario++) {
			long start = ElapsedTime.systemTime();
			double elapsed = 0.0;

			double T = 100;
			double alph = 0.99;
			inputs.setMaxMin();

			//aTest.getLongSim()



			//Creo solucion inicial
			createInitialSolution(scenario);
			System.out.println(initialSolution);

			Solution newSol = new Solution();
			bestSolution = new Solution(initialSolution);
			Solution baseSolution = new Solution(initialSolution);

			//VNS
			while(elapsed < aTest.getMaxTime()){
				int p = P_MIN;
				T = 100;
				while(p <= P_MAX){
					newSol = new Solution(shake(baseSolution,p,scenario));
					if(hasImproved(newSol) || newSol.getTotalCosts()<0){
						System.out.println("I improved shake" + newSol.getTotalScore()+ " "+newSol.getTotalCosts());
					}
					newSol = new Solution(routecache.improve(newSol,inputs,scenario));
					if(hasImproved(newSol) || newSol.getTotalCosts()<0){
						System.out.println("I improved improve" + newSol.getTotalScore()+ " "+newSol.getTotalCosts());
					}
					newSol = new Solution(LocalSearch2(newSol,scenario));
					if(hasImproved(newSol) || newSol.getTotalCosts()<0){
						System.out.println("I improved  LocalSearch2" + newSol.getTotalScore()+ " "+newSol.getTotalCosts());
					}
					newSol = new Solution(LocalSearch4(newSol,scenario));
					if(hasImproved(newSol) || newSol.getTotalCosts()<0){
						System.out.println("I improved LocalSearch4 " + newSol.getTotalScore()+ " "+newSol.getTotalCosts());
					}
					if(hasImproved(newSol)){
						bestSolution = new Solution(newSol);
						baseSolution = new Solution(newSol);
						bestSolution.setTime(elapsed);
						System.out.println("I improved " + newSol.getTotalScore()+ " "+newSol.getTotalCosts());
						p = P_MIN;
					}
					else{ 
						if(SAN(baseSolution,newSol,T)) {
							//System.out.println("Accepted p= "+ p+ " Temp= "+T);
							baseSolution=new Solution(newSol);
							p = P_MIN;
						}
						else {
							p += P_STEP;
						}   
					}
					T=alph*T; 
				}

				elapsed = ElapsedTime.calcElapsed(start, ElapsedTime.systemTime());
			}


			/*** STEP 2*********   */
			bestSolution = new Solution(Stochastic.simulate(bestSolution,aTest.getLongSim(), aTest.getRandomStream(),aTest,inputs,-1,-1));


			cumulativeObjectivefunction[scenario][0]=bestSolution.getTotalScore();
			cumulativeObjectivefunction[scenario][1]=bestSolution.getTotalCosts();
			System.out.println(bestSolution);}
		computingExpectedCostandScore(bestSolution);
		return bestSolution;
	}






	private void computingExpectedCostandScore(Solution bestSolution2) {
		System.out.print(""+ aTest.getLongSim());
		double averageCost=0;
		double averageScore=0;
		double distance=0;
		double score=0;
		for(int scenario=0;scenario<aTest.getLongSim();scenario++) {
			// row=0 <- score // row= 1 <- cost
			score+=cumulativeObjectivefunction[scenario][0];
			distance+=cumulativeObjectivefunction[scenario][1];
		}
		averageCost=distance/aTest.getLongSim();
		averageScore=score/aTest.getLongSim();	
		bestSolution2.setStochCost(averageCost);
		bestSolution2.setStochScore(averageScore);
		System.out.print(" Average cost "+ averageCost+" Average score "+averageScore);
	}





	//Simulated annealing
	private boolean SAN(Solution baseSolution, Solution newSol,double T) 
	{
		double delta=baseSolution.getTotalScore()-newSol.getTotalScore();
		if(Math.abs(delta) > EPSILON) {    		 
			double r = delta / T;
			double a = Math.exp(-r);
			double u = rng.nextDouble();

			if(u < a) return true;            
		}

		return false;
	}



	/**
	 * It creates the initial solution for the VNS metaheuristic, using CWS heuristic, and
	 * determining the best value of alpha to generate the savings list
	 * @param scenario 
	 */
	private void createInitialSolution(int scenario){
		InputsManager.generateDepotEdges(inputs,aTest);
		Solution best = new Solution();
		best.setTotalCosts(Double.MAX_VALUE);


		bestAlpha = 0.0;
		for (double alpha = 0; alpha <= 1; alpha+=0.1){
			InputsManager.generateSavingsList(inputs,aTest,alpha,scenario);
			Solution detSol = new Solution(CWS.solve(inputs,aTest,this.rng,0));
			Collections.sort(detSol.getRoutes());
			if(detSol.getTotalScore() >= best.getTotalScore() ||
					(detSol.getTotalScore() == best.getTotalScore()
					&& detSol.getTotalCosts() < best.getTotalCosts())){
				best = new Solution(detSol);
				bestAlpha = alpha;

			}

		}
		System.out.println("BestAplha: " + bestAlpha);
		initialSolution = best;
	}





	/**
	 * indicates if baseSolution is better than bestSolution
	 * @param baseSolution The solution to compare with BestSolution
	 * @return  True if baseSolution is better than BestSolution, false otherwise
	 */
	private boolean hasImproved(Solution baseSolution){
		double gap =   baseSolution.getTotalCosts() - bestSolution.getTotalCosts();

		return ( (baseSolution.getTotalScore() > bestSolution.getTotalScore()) || 
				((baseSolution.getTotalScore() == bestSolution.getTotalScore()) && (gap  < -0.01)))  ;
	}



	/**
	 * indicates if baseSolution is better than bestSolution
	 * @param baseSolution The solution to compare with BestSolution
	 * @return  True if baseSolution is better than BestSolution, false otherwise
	 */
	private boolean hasImprovedBaseSol(Solution newSol, Solution baseSol){
		double gap =  baseSol.getTotalCosts() - newSol.getTotalCosts();

		return ( (newSol.getTotalScore() > baseSol.getTotalScore()) || 
				((newSol.getTotalScore() == baseSol.getTotalScore()) && (gap > 0.01) )) ;
	}




	private boolean hasImprovedSto(Solution newSol){
		double gap =  bestSolution.getTotalCosts() - newSol.getTotalCosts();

		return (newSol.getStochScore() > bestSolution.getStochScore());
	}




	/**
	 * It destroys a percentage of routes of base Solution and it
	 * generates them again using biased C&W constructive heuristic
	 * @param base  Base solution to shake
	 * @param p percentage of routes of base Solution to delete
	 * @param scenario 
	 */
	private Solution shake(Solution base, int p, int scenario){

		Solution ShakeSol = new Solution(base);
		Solution ShakeSolAux = new Solution(base);

		/*We destroy a percantage of routes*/
		LinkedList<Route> routes = ShakeSol.getRoutes();
		int nRoutesDestroy = (int)  ((routes.size()*p)/100); 

		HashSet<Node> nodesSet = destroyNroutes(routes,nRoutesDestroy);

		/* The nodes destroyed now are not used in the base solution*/
		Set<Node> baseNode = ShakeSol.getNotUsedNodes();
		baseNode.addAll(nodesSet);

		/* We get all nodes not used in the base solution*/
		Node[] nodes = new Node[baseNode.size() + 2];
		nodes[0] = inputs.getNodes()[0];
		nodes[nodes.length - 1] = inputs.getNodes()[inputs.getNodes().length - 1];
		int i = 1;
		for(Node n: baseNode){
			nodes[i] = n;
			++i;
		}


		/* Create a subproblem with previous nodes and solve it with CWS */
		Inputs subInput = new Inputs(inputs,nodes);
		InputsManager.generateDepotEdges(subInput,aTest);

		InputsManager.generateSavingsList(subInput,aTest,bestAlpha,scenario);
		Solution subSol = CWS.solve(subInput,aTest,this.rng,1);
	for(Route r:subSol.getRoutes()) {
		r.updating(inputs, scenario);
		if(r.getCosts()<0) {
			System.out.println("Stop");
		}
	}
		routes.addAll(subSol.getRoutes());

		mergeNotUsedNodes(ShakeSol,subSol);

		/* Merge baseSolution and subSolution and slice routes not used*/

		ShakeSol.sliceSolutionAndSetCost(inputs);

		return ShakeSol;
	}


	/**
	 * It updates the NotUsedNodes set of baseSol using subSol. Every node
	 * used on subSol is removed from the baseSol NotUsedNodes set.
	 * @param baseSol The original solution
	 * @param subSol The sub-solution using nodes not used on baseSol and nodes of the deleted routes
	 */
	private void mergeNotUsedNodes(Solution baseSol, Solution subSol){
		for(Route r: subSol.getRoutes()){
			int i = 0;
			for(Edge e: r.getEdges()){
				if(i != r.getEdges().size())
					baseSol.getNotUsedNodes().remove(e.getEnd());
				++i;
			}
		}
	}

	/**
	 * It destroys n random routes from list routes and it returns a node array
	 * with origin and end positions null and the other positions will be the nodes
	 * that the deleted routes had.
	 * @param routes A list with the total routes
	 * @param n The number of routes to delete
	 * @return the nodes that contained the deleted routes
	 */
	private HashSet<Node> destroyNroutes(List<Route> routes, int n){
		ArrayList<Route> routesDestroy = new ArrayList<>();
		if(!routes.isEmpty()) {
			for (int i = 0; i < n; ++i) {
				int index = this.rng.nextInt(routes.size());
				Route r = routes.get(index);
				routesDestroy.add(r);
			}
		}

		HashSet<Node> nodesD = new HashSet<>();
		for(Route r : routesDestroy){
			int i = 0;
			for(Edge e : r.getEdges()){
				if(i != r.getEdges().size() - 1){
					nodesD.add(e.getEnd());
				}
				++i;
			}
			routes.remove(r);

		}

		return nodesD;
	}




	private  boolean merge(Route route, Edge edge, Node node) {

		boolean ismerge = false;

		Node p = edge.getOrigin();
		Node q = edge.getEnd();

		double cost = (edge.calcCosts(p, node) + edge.calcCosts(node, q) - edge.calcCosts(p, q) );
		if( (route.getCosts() + cost ) <= inputs.gettMax()){
			ismerge = true;
			route.setCosts((route.getCosts() + cost));
		}

		return ismerge;
	}





	/**
	 * It creates the initial solution for the VNS metaheuristic, using CWS heuristic, and
	 * determining the best value of alpha to generate the savings list
	 * @param scenario 
	 */
	private Solution LocalSearch2(Solution Sol, int scenario){
		//Solution LsSolution = new Solution(Sol);
		Boolean isMerge = true;	
		Solution AuxSolLS = new Solution(Sol);
		Solution LsSolution = new Solution(AuxSolLS);
		double min = 0.0;	
		Boolean delRandom = false;


		int nodesToDelete = (int) (inputs.getNodes().length * 0.05);
		boolean stopped = false;
		int ndeleted = 0;
		int thershold = 12;// 
		int del = this.rng.nextInt(3);

		if(del==0){
			min = inputs.getMinprofit();
		}else{
			if(del == 1){
				min = inputs.getMaxprofit()  - thershold;
			}else{
				min = -1.0;
				delRandom = true;
			}

		}

		for(Route r : LsSolution.getRoutes()){
			if(r.getCosts()<0) {
				System.out.println("Stop cost " + r.getCosts());
			}
			
			if(r.getEdges().size() >= 4){
				for(int i = 0;i<r.getEdges().size()-1;i++){
					int index = i;
					if (delRandom == true){
						index = this.rng.nextInt( (int) Math.round(r.getEdges().size() -2));
					}

					Edge edgeIni = r.getEdges().get(index);
					Edge edgeEnd = r.getEdges().get(index+1);

					if( ( (edgeIni.getEnd().getProfit() >= min) && (edgeIni.getEnd().getProfit() <= min + thershold)) || delRandom == true){
						ndeleted++;

						double routeCost = r.getCosts();
						double saveCost = 0;

						int  position = r.getEdges().indexOf(edgeIni); //Posicion del edge n la routa
						r.getEdges().remove(edgeIni); //Borro el egde de la ruta
						saveCost = edgeIni.getCosts();

						int  positionEnd = r.getEdges().indexOf(edgeEnd); //Posicion del edge en la routa
						r.getEdges().remove(edgeEnd); //Borro el egde de la ruta  				
						saveCost +=  edgeEnd.getCosts();
						Edge newEdge =callingEdge(edgeIni.getOrigin(), edgeEnd.getEnd(),scenario);


						r.getEdges().add(position, newEdge);

						double lastScoreRoute = r.getScore();
						if((r.getCosts() - saveCost + newEdge.getCosts())<0) {
							System.out.println("Stop");
						}
						r.setScore(r.getScore() - edgeIni.getEnd().getProfit()); //resto score
						r.setCosts(r.getCosts() - saveCost + newEdge.getCosts());
						if(r.getCosts()<0) {
							System.out.println("Stop");
						}
						
						LsSolution.getNotUsedNodes().add(edgeIni.getEnd());	
						r.setClientsServed(r.getEdges().size()-1); //Customers visitados en la ruta

						LsSolution.setTotalScore(r.getScore() + LsSolution.getTotalScore() - lastScoreRoute); //Actualizo score de la sol
						LsSolution.setTotalCosts(LsSolution.getTotalCosts() - routeCost  + r.getCosts()); //Actualizo cost de la sol (viaje + tiempo de servicio nodo)				

						if(ndeleted == nodesToDelete){
							stopped = true;
						}
					}
				}
				if (stopped == true){
					break;
				}
			}
		}

		AuxSolLS = new Solution(routecache.improve(LsSolution,inputs,scenario)); 
		return AuxSolLS;
	}





	public static Edge callingEdge(Node origin, Node end, int scenario) {
		String key= origin.getId()+","+end.getId();
		Edge newEdge =null; 
		if(inputs.getEdgesDirectory().containsKey(key)) {
			newEdge =inputs.getEdgesDirectory().get(key);
		}
		else {
			newEdge = new Edge(origin, end);
			newEdge.setCosts(newEdge.calcCosts());
			newEdge.setMean(newEdge.getCosts());
			inputs.getEdgesDirectory().put(newEdge.getKey(), newEdge);

		}
		generatingScenarios(newEdge);

		settingCosts(newEdge,scenario);
		return newEdge;
	}





	private static void settingCosts(Edge newEdge, int scenario) {
		if(scenario<2) {
			newEdge.setCosts(newEdge.getMeanCost());
		}
		else {
			newEdge.setCosts(newEdge.getStoCosts(scenario));
		}
	}





	private static void generatingScenarios(Edge e) {
		if(e.getStoCost()==null) {
			double[] times= new double[aTest.getLongSim()];
			if(e.getMeanCost()==0) {
				System.out.print("Stop");
			}
			double mean = e.getMeanCost(); // esta es la distancia de cada eje 
			double newArc=mean;
			for(int i=0;i<Math.max(aTest.getLongSim(),1);i++) {
				if(mean>0) {
					newArc =getStochasticValue(aTest.getRandomStream(),mean,aTest.getVariance());	
					times[i]=newArc;
				}
			}	
			e.setScenarios(times);
		}	

	}





	public static double getStochasticValue(RandomStream stream, double mean, float variance) {
		double squareSigma = Math.log(1 + (variance / Math.pow(mean, 2)));
		double mu = Math.log(mean) - squareSigma / 2;
		double sigma = Math.sqrt(squareSigma);
		return LognormalGen.nextDouble(stream, mu, sigma);
	}






	/**
	 * It creates the initial solution for the VNS metaheuristic, using CWS heuristic, and
	 * determining the best value of alpha to generate the savings list
	 * @param scenario 
	 */
	private Solution LocalSearch4(Solution Sol, int scenario){
		//Solution LsSolution = new Solution(Sol);
		Boolean isMerge = true;	
		Solution AuxSolLS = new Solution(Sol);
		Solution LsSolution = new Solution(AuxSolLS);
		Solution BaseSol = new Solution(AuxSolLS);
		Boolean improve = true;

		while(improve == true){
			for(Route r : LsSolution.getRoutes()){
				for(int i = 0;i<r.getEdges().size();i++){

					int index = i;
					Edge edge = r.getEdges().get(index);

					//calculate function evaluation for all the availables nodes for the selected edge
					//caculate distance and order between nodes(j's) and Edge(p,q)  return list in increasing order
					if(LsSolution.getNotUsedNodes().size() == 0){
						isMerge = false;
						break; // there aren't free nodes (all in use)
					}
					TreeSet<PairBest> vectDist = MultiStart.calEvaluationFunction(LsSolution.getNotUsedNodes(),edge,r,inputs);

					index = MultiStart.getRandomPosition(0.30, rng, vectDist.size()); //select node using Bias randomization

					PairBest b = MultiStart.obtainDistanceNode(index,vectDist);	 

					if(merge(r,edge,b.getvalue())){

						int  position = r.getEdges().indexOf(edge); //Posicion del edge en la routa
						r.getEdges().remove(edge); //Borro el egde de la ruta
						//System.out.println("Mejoraaaa LS");

						//P -> J
						Edge pjEdge = callingEdge(edge.getOrigin(), b.getvalue(), scenario);

						//pjEdge.setCosts(pjEdge.calcCosts(edge.getOrigin(), b.getvalue()));
						r.getEdges().add(position, pjEdge);

						// System.out.println(b.getvalue().getProfit() + " " + route.getScore());
						double lastScoreRoute = r.getScore();
						r.setScore(b.getvalue().getProfit() + r.getScore()); //Sumo score

						//J -> Q
						Edge jqEdge =callingEdge(b.getvalue(), edge.getEnd(), scenario);
						r.getEdges().add(position + 1, jqEdge);

						r.setClientsServed(r.getEdges().size()-1); //Customers visitados en la ruta
						if(r.getCosts()<0) {
							System.out.println("Stop");
						}
						
						r.updating(inputs,scenario);
						LsSolution.setTotalScore(r.getScore() + LsSolution.getTotalScore() - lastScoreRoute); //Actualizo score de la sol
						LsSolution.setTotalCosts(LsSolution.getTotalCosts() + r.getCosts()); //Actualizo cost de la sol (viaje + tiempo de servicio nodo)
						if(LsSolution.getTotalCosts()<0) {
							System.out.println("Stop");
						}
						LsSolution.getNotUsedNodes().remove(b.getvalue());	

						if( (LsSolution.getTotalScore() > AuxSolLS.getTotalScore()) ||
								( (LsSolution.getTotalScore() == AuxSolLS.getTotalScore()) && (LsSolution.getTotalCosts() < AuxSolLS.getTotalCosts()))){
							AuxSolLS = new Solution(routecache.improve(LsSolution,inputs,scenario)); 
						}
					}	
				} 
			}

			if ((AuxSolLS.getTotalScore() > BaseSol.getTotalScore()) ||
					( (AuxSolLS.getTotalScore() == BaseSol.getTotalScore()) && (AuxSolLS.getTotalCosts() < BaseSol.getTotalCosts()))){
				LsSolution = new Solution(routecache.improve(AuxSolLS,inputs,scenario)); 
				BaseSol = new Solution(routecache.improve(LsSolution,inputs,scenario)); 
			}else{improve = false;}
		}

		
		AuxSolLS.sliceSolutionAndSetCost(inputs);
		if(hasImproved(AuxSolLS)){
			System.out.println("I improved LocalSearch4 " + AuxSolLS.getTotalScore()+ " "+AuxSolLS.getTotalCosts());
		}
		return AuxSolLS;
	}



}
