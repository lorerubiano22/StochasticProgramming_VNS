package alg;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Locale;

import ilog.concert.IloException;
import ilog.concert.IloLinearNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class TOP_deterministicVersion1 {
	public static void solveMe(Inputs inp, Test t) {

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/*SETS*/
		int V=inp.getVehNumber(); // set of vehicles
		int K=t.getLongSim()+1; // set of scenarios
		int n=inp.getNodes().length; //sset of node
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/*PARAMETERS*/


		//Traveling time
		double [][][]TV= new double [n][n][K];
		for(Edge e:inp.getedgeList()) {
			int indexI=e.getOrigin().getId();
			int indexJ=e.getEnd().getId();
			for(int k=0;k<e.getStchTV().size();k++) {
				double tv=e.getStchTV().get(k);
				TV[indexI][indexJ][k]=tv;
			}
		}

		//Rewards
		double []U= new double [n];
		for(int i=0; i<n;i++) {
			Node IDnode=inp.getNodes()[i];
			U[i]=IDnode.getProfit();}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/*MODEL*/
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		try {
			IloCplex cplex= new IloCplex();

			/*DECISION VARIABLES*/

			//Binary variable X(i,j,d) if UAV travel from i to j using vehicle d
			IloNumVar[][][] x= new IloNumVar[n][n][V];
			for(int i=0;i<n;i++){
				for(int j=0;j<n;j++){ 
					x[i][j]=cplex.boolVarArray(V);
				}
			}

			//Auxiliary Variable to eliminate subtours. it denotes the position of the node i in the route assigned to vehicle d
			IloNumVar[][] w=new IloNumVar[n][V];
			for(int i=0;i<n;i++){
				w[i]=cplex.numVarArray(V, 0,Double.MAX_VALUE);
			}

			// outcome variables to route length
			//Length of routes d in the scenario k 
			IloNumVar[][] lenght=new IloNumVar[V][K];
			for(int d=0;d<V;d++){
				lenght[d]=cplex.numVarArray(K, 0,Double.MAX_VALUE);
			}

			/*OBJECTIVE FUNCTION*/

			IloLinearNumExpr obj = cplex.linearNumExpr();
			for(int i=0; i<n;i++) {
				for(int j=0;j<n;j++ ) {
					if(i!=j) {
						for(int d=0;d<inp.getVehNumber();d++ ) {
							obj.addTerm(U[i], x[i][j][d]);
						}
					}
				}
			}
			cplex.addMaximize(obj);

			/*CONSTRAINTS*/
			
			
			//1. Todos los drones parten del depot inicial origen_cust(d)..sum((j)$(ord(j)>=1 and ord(j)<end),X('0',j,d))=l=1;  
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				for(int j= 1; j<n;j++) {
					expr2.addTerm(1.0, x[0][j][d]);
				}
				cplex.addEq(expr2,1.0);
			}
			
			
			

			//2. Todos los drones tienen que llegar al depot final cust_depot(d)..sum((i,j)$(ord(i)>=0 and ord(i)<end-1 and ord(j)=end),X(i,n,d))=l=1;
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr expr = cplex.linearNumExpr();
				for(int i= 0; i<n-1;i++) {
					expr.addTerm(1.0, x[i][(n-1)][d]);
				}
				cplex.addEq(expr,1.0);
			}



		

//			//3. no existe conección entre los depots origen_end(j,d)$(ord(j)=end)..,X('0',n,d)=l=1; 
//			for(int d=0;d<inp.getVehNumber();d++) {
//				IloLinearNumExpr expr2 = cplex.linearNumExpr();
//				expr2.addTerm(1.0, x[0][(n-1)][d]);
//				cplex.addEq(expr2,0.0);
//			}
			
			

			//4. cada punto j tiene un solo origen i visit_cust(j)$(ord(j)<end and ord(j)>1)..sum((i,d),X(i,j,d))=l=1; 
			for(int j=1;j<n;j++) {
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				for(int i=0;i<n-1;i++) {
					for(int d=0;d<V;d++) {
						expr2.addTerm(1.0, x[i][j][d]);}}
				cplex.addLe(expr2,1);
			}
			
			//4. cada punto i tiene un solo destino j visit_cust(j)$(ord(j)<end and ord(j)>1)..sum((i,d),X(i,j,d))=l=1; 
			for(int i=0;i<n-1;i++) {
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				for(int j=1;j<n;j++) {
					for(int d=0;d<V;d++) {
						expr2.addTerm(1.0, x[i][j][d]);}}
				cplex.addLe(expr2,1);
			}


			
			//3. flow(o,d)$(ord(o)<>1 and ord(o)<>end)..sum(i,X(i,o,d))=e=sum(j,X(o,j,d)); 
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int o=1; o<(n-1);o++) {
					IloLinearNumExpr expr3 = cplex.linearNumExpr();
					for(int i= 0; i<n-1;i++) {
						for(int j= 1; j<n;j++) {
							expr3.addTerm(1.0, x[i][o][d]);
							expr3.addTerm(-1.0, x[o][j][d]);
						}
					}
					cplex.addEq(expr3, 0); 
				}

			}


//
//			// 8. Subtour(i,j,d)$(ord(i)>1 and ord(i)<end and ord(j)>1 and ord(j)<end).. U(i,d)-U(j,d)+(card(i)*X(i,j,d))=l=card(i)-1; 
//			for(int d=0;d<inp.getVehNumber();d++) {
//				for(int i=0;i<n;i++) {
//					for(int j= 0; j<n;j++) {
//						IloLinearNumExpr expr8 = cplex.linearNumExpr();
//						expr8.addTerm(1.0, w[i][d]);
//						expr8.addTerm(-1.0, w[j][d]);
//						expr8.addTerm(n, x[i][j][d]);
//						cplex.addLe(expr8, n-1);
//					}
//				}
//			}
			

			//9. LenghtRoute(d,k)..Lenght(d,k)=e=sum((i,j),X(i,j,d)*(TV(i,j)));
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int k=0;k<K;k++) {
					IloLinearNumExpr expr10 = cplex.linearNumExpr();
					for(int i= 0; i<n;i++) {
						for(int j= 0; j<n;j++) {
							expr10.addTerm(TV[i][j][k], x[i][j][d]);
						}
					}
					cplex.addLe(expr10,inp.gettMax());
				}
			}

			
			
			/*SOLVING*/
			cplex.setParam(IloCplex.Param.TimeLimit, 10000);
			if (cplex.solve()) {
				String CPLEX_file= new String(t.getInstanceName()+t.getLongSim()+"_"+"_"+t.getVariance()+"Stoch_CPLEX.txt");
				writeFile(CPLEX_file,cplex,x,V,n,K);
				//PrintWriter bw = new PrintWriter(CPLEX_file);
				System.out.println("obj = " + cplex.getObjValue());
				for (int d=0; d<V; d++){
					System.out.println("Route "+d+1);
					for (int i=0; i<n; i++){
						for (int j=0; j<n; j++){
							if (i!=j && cplex.getValue(x[i][j][d])>0.9999) {
								System.out.println(i+" - "+j);
							}
						}
					}
					System.out.println("\n");
				}


			}
			else {
				System.out.println("Model not solved");
			};
			cplex.end(); 

		}
		catch(IloException e) {
			e.printStackTrace();
		}
	}
	static void writeFile(String filename, IloCplex cplex, IloNumVar[][][] x,int V, int n, int K) {
		try {
			PrintWriter bw = new PrintWriter(filename);
			try {
				bw.printf(Locale.US, "obj = "+ cplex.getObjValue());
				bw.println();
				for (int d=0; d<V; d++){
					bw.println("Route "+d+1);
					for (int i=0; i<n; i++){
						for (int j=0; j<n; j++){
							if (i!=j && cplex.getValue(x[i][j][d])>0.9999) {
								bw.println(i+" - "+j);
							}
						}
					}
					bw.println();
					bw.print(cplex.getCplexTime());
				}
			}
			catch(IloException e) {
				e.printStackTrace();
			}
			bw.println();
			bw.flush();
		} 
		catch (IOException e) {
			//why does the catch need its own curly?
		}
	}
}