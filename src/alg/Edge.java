package alg;
import java.io.Serializable;
import java.util.Comparator;

import umontreal.iro.lecuyer.randvar.LognormalGen;
import umontreal.iro.lecuyer.rng.RandomStream;

/**
 * @author Angel A. Juan - ajuanp(@)gmail.com
 * @version 130112
 */
public class Edge implements Serializable
{
	/* INSTANCE FIELDS & CONSTRUCTOR */
	private String key=""; // origin node
	private Node origin; // origin node
	private Node end; // end node
	private double meanCost = 0.0; // edge costs
	private double cost = 0.0; // cost of the edge accordoing the scenario
	private double savings = 0.0; // edge savings (Clarke & Wright)
	private double classicSavings = 0.0;
	private Route inRoute = null; // route containing this edge (0 if no route assigned)
	private Edge inverseEdge = null; // edge with inverse direction
	private double[] stoCosts =null;; // edge costs




	public Edge(Node originNode, Node endNode, int scenarios) 
	{   origin = originNode;
	end = endNode;
	stoCosts= new double[scenarios];
	key=origin.getId()+","+endNode.getId();
	}

	public Edge(Edge e){   
		this.key=e.getKeyEdge();
		this.origin = e.origin;
		this.end = e.end; 
		this.cost = e.cost;
		this.meanCost = e.meanCost;
		this.savings = e.savings; 
		this.classicSavings = e.classicSavings;
		if(e.inRoute !=null){
			this.inRoute = new Route (e.inRoute);
		}else{
			this.inRoute = null;	
		}
	}



	/* SET METHODS */
	public void setMeanCosts(double c){meanCost = c;}
	public void setCosts(double c){cost = c;}
	public void setSavings(double s){savings = s;}
	public void setInRoute(Route r){inRoute = r;}
	public void setInverse(Edge e){inverseEdge = e;}

	/* GET METHODS */
	public Node getOrigin(){return origin;}
	public Node getEnd(){return end;}
	public double getMeanCosts(){return meanCost;}
	public double getCosts(){return cost;}
	public double getSavings(){return savings;}
	public Route getInRoute(){return inRoute;}
	public Edge getInverseEdge(){return inverseEdge;}
	public double getClassicSavings() {return classicSavings;}
	public double getTimeScenario(int i) {return stoCosts[i];}
	public String getKeyEdge() {return key;}


	/* AUXILIARY METHODS */

	public double calcCosts()
	{   double X1 = origin.getX();
	double Y1 = origin.getY();
	double X2 = end.getX();
	double Y2 = end.getY();
	double d = Math.sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1));
	return d;
	}

	public static double calcCosts(Node origin, Node end)
	{   double X1 = origin.getX();
	double Y1 = origin.getY();
	double X2 = end.getX();
	double Y2 = end.getY();
	double d = Math.sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1));
	return d;
	}

	public double calcSavings(Node origin, Node end)
	{
		// Costs of origin depot to end node
		double Coj = calcCosts(origin,this.end);
		// Costs of origin node to end depot
		double Cie = calcCosts(this.origin,end);
		// Costs of originNode to endNode
		double Cij = meanCost;

		//Return cost depot to savings
		return Coj + Cie - Cij;
	}


	public double calcSavings(Node origin, Node end, double alpha)
	{
		// Costs of origin depot to end node
		double Coj = calcCosts(origin,this.end);
		// Costs of origin node to end depot
		double Cie = calcCosts(this.origin,end);
		// Costs of originNode to endNode
		double Cij = meanCost;

		//Return cost depot to savings
		double Sij = Coj + Cie - Cij;
		classicSavings = Sij;
		return alpha*Sij + (1-alpha)*(this.origin.getProfit() + this.end.getProfit());
	}

	static final Comparator<Edge> minDistance = new Comparator<Edge>(){
		public int compare(Edge a1, Edge a2){
			if (a1.meanCost < a2.meanCost) return 1;
			if (a1.meanCost > a2.meanCost) return -1;
			return 0;
		}
	};

	static final Comparator<Edge> savingsComp = new Comparator<Edge>(){
		public int compare(Edge a1, Edge a2){
			if (a1.savings < a2.savings) return 1;
			if (a1.savings > a2.savings) return -1;
			return 0;
		}
	};





	@Override
	public String toString() 
	{ 
		String s = "";
		s = s.concat("\nEdge origin: " + this.getOrigin());
		s = s.concat("\nEdge end: " + this.getEnd());
		s = s.concat("\nEdge costs: " + (this.getCosts()));
		s = s.concat("\nEdge savings: " + (this.getSavings()));
		return s;
	}

	public void setTimeByScenario(Test aTest) {
		double mean = this.getMeanCosts(); // esta es la distancia de cada eje 
		double newArc=mean;
		for(int i=0;i<aTest.getLongSim();i++) {
			if(mean>0) {
				newArc =getStochasticValue(aTest.getRandomStream(),mean,aTest.getVariance());	
				stoCosts[i]=newArc;
			}
		}
	}

	public double getStochasticValue(RandomStream stream, double mean, float variance) {
		double squareSigma = Math.log(1 + (variance / Math.pow(mean, 2)));
		double mu = Math.log(mean) - squareSigma / 2;
		double sigma = Math.sqrt(squareSigma);
		return LognormalGen.nextDouble(stream, mu, sigma);

	}

}