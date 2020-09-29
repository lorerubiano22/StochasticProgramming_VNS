package alg;

import java.util.*;

/**
 * Created by lluc on 2/06/17.
 * This class encapsulates the metaheuristic VNS
 */
public class VNS {

	private static final int P_MIN = 1 ;
	private static final int P_MAX = 100 ;
	private static final int P_STEP = 2 ;
	private Test aTest;
	private Inputs inputs;
	private Random rng;
	private double bestAlpha;
	private Solution initialSolution;
	private Solution bestSolution;
	private static final double EPSILON = 10E-6;
	RouteCache  routecache;
	private double[][] cumulativeObjectivefunction=null;


	VNS(Test myTest, Inputs myInputs, Random myrng){
		aTest = myTest;
		inputs = myInputs;
		rng = myrng;
		cumulativeObjectivefunction=new double[aTest.getLongSim()][2];
		routecache = new RouteCache(inputs);
	}

	/**
	 * It applies VNS metaheuristic to solve the problem specified at attribute inputs, using test instance class aTest
	 * @return The best solution found by VNS
	 */
	public Solution solveDet(){
		for(int scenario=0;scenario<aTest.getLongSim();scenario++  ) { // aca es donde se va a iterar sobre cada simulción
			//	for(int scenario=0;scenario<aTest.getLongSim();scenario++  ) { // aca es donde se va a iterar sobre cada simulción
			//	int scenario=0; // to remove
			long start = ElapsedTime.systemTime();
			double elapsed = 0.0;
			double T = 100;
			double alph = 0.99;
			inputs.setMaxMin();
			//Creo solucion inicial
			createInitialSolution(scenario); //
			System.out.println(initialSolution);
			Solution newSol = new Solution();
			bestSolution = new Solution(initialSolution);
			Solution baseSolution = new Solution(initialSolution);
			//VNS
			while(elapsed < aTest.getMaxTime()){
				int p = P_MIN;
				T = 100;
				while(p <= P_MAX){
					if(p==51 && elapsed==22.3721493) {
						System.out.println("Stop");
					}
					newSol = new Solution(shake(baseSolution,p,scenario ));
					if(newSol.getTotalCosts()<0) {
						System.out.println("To check new Sol" +"iteration" + p + bestSolution.getTotalScore()+ " "+newSol.getTotalCosts());

					}
					newSol = new Solution(routecache.improve(newSol));
					//checkingDistances(newSol);
					if(newSol.getTotalCosts()<0) {
						System.out.println("To check new Sol" +"iteration" + p + bestSolution.getTotalScore()+ " "+newSol.getTotalCosts());

					}
					newSol = new Solution(LocalSearch2(newSol));
					//	checkingDistances(newSol);
					if(newSol.getTotalCosts()<0) {
						System.out.println("To check new Sol" +"iteration" + p + bestSolution.getTotalScore()+ " "+newSol.getTotalCosts());

					}
					newSol = new Solution(LocalSearch4(newSol));
					//		checkingDistances(newSol);
					if(newSol.getTotalCosts()<0) {
						System.out.println("To check new Sol" +"iteration" + p + bestSolution.getTotalScore()+ " "+newSol.getTotalCosts());

					}
					if(hasImproved(newSol)){
						bestSolution = new Solution(newSol);
						if(bestSolution.getTotalCosts()<0) {
							System.out.println("To check bestSolution" +"iteration" + p + bestSolution.getTotalScore()+ " "+newSol.getTotalCosts());
						}
						baseSolution = new Solution(newSol);

						bestSolution.setTime(elapsed);
						System.out.println("I improved " + bestSolution.getTotalScore()+ " "+bestSolution.getTotalCosts());
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
			if(bestSolution.getTotalCosts()<0) {
				System.out.println("To check bestSolution after simulation"+ bestSolution.getTotalScore()+ " "+newSol.getTotalCosts());
			}
			System.out.println(bestSolution);
			cumulativeObjectivefunction[scenario][0]=bestSolution.getTotalScore();
			cumulativeObjectivefunction[scenario][1]=bestSolution.getTotalCosts();
		}
		computeExpectedScoreBySAA();
		return bestSolution;
	}






	private void computeExpectedScoreBySAA() {

		// [row=0] es el score
		// [row=2] es el 
		double scenarioScore=0;
		double scenarioDistance=0;

		for(int i=0;i<aTest.getLongSim();i++) {
			scenarioDistance+=cumulativeObjectivefunction[i][1];
			scenarioScore+=cumulativeObjectivefunction[i][0];
		}		
		double expectedScore=scenarioScore/aTest.getLongSim();
		double expectedDistance=scenarioDistance/aTest.getLongSim();
		bestSolution.setScoreSAA(expectedScore);
		bestSolution.setDistanceSAA(expectedDistance);

	}

	/**
	 * It applies SimVNS metaheuristic to solve the problem specified at attribute inputs, using test instance class aTest
	 * @return The best solution found by VNS
	 */






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
		Solution best = new Solution();
		best.setTotalCosts(Double.MAX_VALUE);
		// 1. setting los costos para cada de acuerdo al escenario
		changingTimesAccorgingScenario(inputs,scenario);
		for(Edge e:inputs.getEdgesDirectory().values()) {
			System.out.println("Edge_"+e.getKeyEdge()+"_mean_"+e.getMeanCosts()+"_timeScenario_1___"+e.getCosts()+"_timeScenario__3_"+e.getTimeScenario(3));
		}
		bestAlpha = 0.0;
		for (double alpha = 0; alpha <= 1; alpha+=0.1){
			InputsManager.generateSavingsList(inputs,alpha,aTest, scenario);				
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
		inputs.getSavingsDirectory().put(scenario, inputs.getSavings());	
		initialSolution = best;
	}


	private void changingTimesAccorgingScenario(Inputs inputs, int scenario) {
		for(Edge e:inputs.getEdgesDirectory().values()) {
			e.setCosts(e.getTimeScenario(scenario));
		}
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
		generateDepotEdges(subInput,aTest);

		//generateSavingsList(Inputs inputs,double alpha, Test t, int scenario)
		generateSavingsList(subInput,bestAlpha,aTest,scenario);
		Solution subSol = CWS.solve(subInput,aTest,this.rng,1);
		routes.addAll(subSol.getRoutes());
		mergeNotUsedNodes(ShakeSol,subSol);
		//	checkingDistances(subSol);
		/* Merge baseSolution and subSolution and slice routes not used*/

		ShakeSol.sliceSolutionAndSetCost(inputs);
		//	checkingDistances(ShakeSol);
		System.out.println("shake");
		return ShakeSol;

	}


	private void checkingDistances(Solution shakeSol) {
		// 1. Llamar las rutas de la solución 
		int route=-1;
		double computedSolCost=0;
		for(Route r:shakeSol.getRoutes()) {
			route++;
			// 2. Llamar los ejes de cada ruta
			double distRoute=0;
			for(Edge e:r.getEdges()) {
				// 3. Calcular la distancia de cada ruta
				distRoute+=e.getCosts();
			}
			// 4. Verificar que todo coincida
			System.out.println("Route_"+route+" distance compute "+distRoute+" current route distance "+ r.getCosts());		
			// 5. solution cost
			computedSolCost+=distRoute;
			if(Math.round(distRoute)!=Math.round(r.getCosts())) {
				System.out.println("Route Cost_wrong");
			}	
		}

		System.out.println("Solution Cost_"+ "computed cost "+computedSolCost+" current route distance "+ shakeSol.getTotalCosts());	
		if(computedSolCost!=shakeSol.getTotalCosts()) {
			System.out.println("Solution Cost_wrong");
		}	
	}

	private void generateSavingsList(Inputs subInput, double alpha, Test aTest2, int scenario) {
		// 2. compute los savings para cada eje
		int nNodes = subInput.getNodes().length;
		LinkedList<Edge> savingsList = new LinkedList<Edge>();
		Node origin = subInput.getNodes()[0];
		Node end = subInput.getNodes()[nNodes - 1];
		for( int i = 1; i < nNodes - 1; i++ ) {// node 0 is the depot
			for( int j = i + 1; j < nNodes - 1; j++ ){
				Node iNode = subInput.getNodes()[i];
				Node jNode = subInput.getNodes()[j];
				// computing the saving for edge (i)---- (j)
				// Costs of origin depot to end node
				Edge ijEdge =inputs.getEdge(iNode, jNode);
				if(ijEdge==null) {
					ijEdge = new Edge(iNode, jNode,aTest.getLongSim());
					ijEdge.setMeanCosts(ijEdge.calcCosts());
					ijEdge.setCosts(ijEdge.calcCosts());
					ijEdge.setTimeByScenario(aTest);
					subInput.getEdgesDirectory().put(ijEdge.getKeyEdge(), ijEdge);
				}
				// Costs of origin depot to end node
				double Coj = inputs.getEdge(origin, jNode).getTimeScenario(scenario);// se calculan los savins de acuerdo con el escenarios
				// Costs of origin node to end depot
				double Cie = inputs.getEdge(iNode, end).getTimeScenario(scenario);
				// Costs of originNode to endNode
				double Cij = ijEdge.getTimeScenario(scenario);
				// computing savings
				double classicSaving = Coj + Cie - Cij;
				double weightedSaving = alpha*classicSaving + (1-alpha)*(iNode.getProfit() + jNode.getProfit());
				ijEdge.setSavings(weightedSaving);
				// inverse edge
				Edge jiEdge =inputs.getEdge(jNode, iNode);
				if(jiEdge==null) {
					jiEdge = new Edge(jNode, iNode,aTest.getLongSim());
					jiEdge.setMeanCosts(jiEdge.calcCosts());
					jiEdge.setCosts(jiEdge.calcCosts());
					jiEdge.setTimeByScenario(aTest);
					subInput.getEdgesDirectory().put(jiEdge.getKeyEdge(), jiEdge);
				}
				// Costs of origin depot to end node
				Coj = inputs.getEdge(origin, iNode).getTimeScenario(scenario);// se calculan los savins de acuerdo con el escenarios
				// Costs of origin node to end depot
				Cie = inputs.getEdge(jNode, end).getTimeScenario(scenario);
				// Costs of originNode to endNode
				Cij = jiEdge.getTimeScenario(scenario);
				// computing savings
				classicSaving = Coj + Cie - Cij;
				weightedSaving = alpha*classicSaving + (1-alpha)*(iNode.getProfit() + jNode.getProfit());
				jiEdge.setSavings(weightedSaving);	
				savingsList.add(ijEdge);
				savingsList.add(jiEdge);
			}
		}
		// Construct the savingsList by sorting the edgesList. Uses the compareTo()
		//  method of the Edge class (TIE ISSUE #1).
		savingsList.sort(Edge.savingsComp);
		subInput.setList(savingsList);

	}

	//inputs.getEdge(
	public  void generateDepotEdges(Inputs subInp, Test aTest){
		Node[] nodes = inputs.getNodes();
		Node depot = nodes[0]; // depot is always node 0
		Node end = nodes[nodes.length-1];
		// Create diEdge and idEdge, and set the corresponding costs
		int nNodes = subInp.getNodes().length;
		LinkedList<Edge> distanceList = new LinkedList<Edge>();
		for( int i = 1; i < nodes.length - 1; i++ ) // node 0 is depot
		{   
			Node iNode = nodes[i];
			Edge diEdge = inputs.getEdge(depot, iNode);
			distanceList.add(diEdge);
			//Edge diEdge = new Edge(depot, iNode, aTest.getLongSim());
			//
			//Edge idEdge = new Edge(iNode, end,aTest.getLongSim());  // solo se crea el eje y el array con los parametros
			Edge idEdge = inputs.getEdge(iNode, end);

			distanceList.add(idEdge);
		}
		distanceList.sort(Edge.minDistance);
		subInp.setdistanceDepot(distanceList);      
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
	 */
	private Solution LocalSearch2(Solution Sol){
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
						Edge newEdge =inputs.getEdge(edgeIni.getOrigin(), edgeEnd.getEnd());
						//Edge newEdge = new Edge(edgeIni.getOrigin(), edgeEnd.getEnd());

						//newEdge.setCosts(newEdge.calcCosts(edgeIni.getOrigin(), edgeEnd.getEnd()));
						r.getEdges().add(position, newEdge);

						double lastScoreRoute = r.getScore();
						r.setScore(r.getScore() - edgeIni.getEnd().getProfit()); //resto score
						r.setCosts(r.getCosts() - saveCost + newEdge.getCosts());

						LsSolution.getNotUsedNodes().add(edgeIni.getEnd());	
						r.setClientsServed(r.getEdges().size()-1); //Customers visitados en la ruta

						//LsSolution.setTotalScore(r.getScore() + LsSolution.getTotalScore() - lastScoreRoute); //Actualizo score de la sol
						//LsSolution.setTotalCosts(LsSolution.getTotalCosts() - routeCost  + r.getCosts()); //Actualizo cost de la sol (viaje + tiempo de servicio nodo)				
						LsSolution.sliceSolutionAndSetCost(inputs);
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

		AuxSolLS = new Solution(routecache.improve(LsSolution)); 
		System.out.println("LocalSearch2");
		return AuxSolLS;
	}





	/**
	 * It creates the initial solution for the VNS metaheuristic, using CWS heuristic, and
	 * determining the best value of alpha to generate the savings list
	 */
	private Solution LocalSearch4(Solution Sol){
		//Solution LsSolution = new Solution(Sol);
		Boolean isMerge = true;	
		Solution AuxSolLS = new Solution(Sol);
		Solution LsSolution = new Solution(AuxSolLS);
		Solution BaseSol = new Solution(AuxSolLS);
		Boolean improve = true;
		//checkingDistances(LsSolution);
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
						// información de la ruta
						double lastScoreRoute = r.getScore();
						double lastDistance=r.getCosts();
						// información de la solución
						double lastScoreSol = LsSolution.getTotalScore()-lastScoreRoute;
						double lastDistanceSol=LsSolution.getTotalCosts()-lastDistance;
						// setting information of solution
						LsSolution.setTotalScore(lastScoreSol); //Actualizo score de la sol
						LsSolution.setTotalCosts(lastDistanceSol); //Actualizo cost de la sol (viaje + tiempo de servicio nodo)




						// route
						int  position = r.getEdges().indexOf(edge); //Posicion del edge en la routa
						r.getEdges().remove(edge); //Borro el egde de la ruta
						//System.out.println("Mejoraaaa LS");

						//P -> J
						Edge pjEdge = inputs.getEdge(edge.getOrigin(), b.getvalue());
						r.getEdges().add(position, pjEdge);

						//System.out.println("Route_"+ r.toString());	
						//	r.addCosts(pjEdge);
						// System.out.println(b.getvalue().getProfit() + " " + route.getScore());

						r.setScore(b.getvalue().getProfit() + r.getScore()); //Sumo score

						//J -> Q
						Edge jqEdge = inputs.getEdge(b.getvalue(), edge.getEnd());
						r.getEdges().add(position + 1, jqEdge);
						r.setClientsServed(r.getEdges().size()-1); //Customers visitados en la ruta
						//System.out.println("Route_"+ r.toString());	


						// computing Distances Route
						computeDistanceScoreRoute(r);
						LsSolution.sliceSolutionAndSetCost(inputs);

						LsSolution.getNotUsedNodes().remove(b.getvalue());	
						//checkingDistances(LsSolution);
						if( (LsSolution.getTotalScore() > AuxSolLS.getTotalScore()) ||
								( (LsSolution.getTotalScore() == AuxSolLS.getTotalScore()) && (LsSolution.getTotalCosts() < AuxSolLS.getTotalCosts()))){
							AuxSolLS = new Solution(routecache.improve(LsSolution)); 
						}
					}	
				} 
			}

			if ((AuxSolLS.getTotalScore() > BaseSol.getTotalScore()) ||
					( (AuxSolLS.getTotalScore() == BaseSol.getTotalScore()) && (AuxSolLS.getTotalCosts() < BaseSol.getTotalCosts()))){
				LsSolution = new Solution(routecache.improve(AuxSolLS)); 
				BaseSol = new Solution(routecache.improve(LsSolution)); 
			}else{improve = false;}
		}
		return AuxSolLS;
	}

	private void computeDistanceScoreRoute(Route r) {
		double distRoute=0;
		double score=r.getEdges().get(0).getOrigin().getProfit();
		for(Edge e:r.getEdges()) {
			// 3. Calcular la distancia de cada ruta
			distRoute+=e.getCosts();
			score+=e.getEnd().getProfit();
		}
		// 4. Verificar que todo coincida

		System.out.println("LocalSearch4");		
		r.setScore(score);
		r.setCosts(distRoute);		
		System.out.println("Route_"+" distance compute "+distRoute+" current route distance "+ r.getCosts());	
	}



}
