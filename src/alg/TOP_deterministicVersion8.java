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

public class TOP_deterministicVersion8 {
	public static void solveMe(Inputs inp, Test t) {

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/*SETS*/
		int V=inp.getVehNumber(); // set of vehicles
		int K=t.getLongSim(); // set of scenarios
		int n=inp.getNodes().length; //sset of node
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/*PARAMETERS*/



		System.out.println("Num veh:" +inp.getVehNumber());
		System.out.println("Dr:" + inp.gettMax());
		System.out.println("K:" + K);

		//Traveling time
		double [][]TV= new double [n][n];
		for(Edge e:inp.getedgeList()) {
			int indexI=e.getOrigin().getId();
			int indexJ=e.getEnd().getId();

			double tv=e.getCosts();
			TV[indexI][indexJ]=tv;

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
			IloNumVar[][][] x= new IloNumVar[n][n][V+1];
			for(int i=0;i<n;i++){
				for(int j=0;j<n;j++){ 
					x[i][j]=cplex.boolVarArray(V+1);
				}
			}

			//Auxiliary Variable to eliminate subtours. it denotes the position of the node i in the route assigned to vehicle d
			IloNumVar[][] w=new IloNumVar[n][V+1];
			for(int i=0;i<n;i++){
				w[i]=cplex.numVarArray(V, 0,Double.MAX_VALUE);
			}



			// outcome variables to route length
			//Length of routes d 
			IloNumVar[] lenght=new IloNumVar  [inp.getVehNumber()];
			for(int d=0;d<inp.getVehNumber();d++){
				lenght[d]=cplex.numVar(0, Double.MAX_VALUE);
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
			/********************/
			/********************/

			/*A. DEPOT*/


			/********************/
			//1. El número de drones que salen del depot no puede exceder la cant disponible
			IloLinearNumExpr C1 = cplex.linearNumExpr();
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int j= 1; j<n-1;j++) {
					C1.addTerm(1.0, x[0][j][d]);
					cplex.addLe(C1, inp.getVehNumber());
				}
			}

			IloLinearNumExpr C1a = cplex.linearNumExpr();
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int i= 1; i<n-1;i++) {
					C1a.addTerm(1.0, x[i][n-1][d]);
					cplex.addLe(C1a, inp.getVehNumber());
				}
			}
			/********************/


			/********************/
			//2. Las rutas empiezan y terminan en el depot
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr C2 = cplex.linearNumExpr();
				for(int j= 0; j<n;j++) { 
					C2.addTerm(1.0, x[0][j][d]);
				}
				cplex.addEq(C2,1.0);
			}

			//terminan
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr c21 = cplex.linearNumExpr();
				for(int i= 0; i<n;i++) {
					c21.addTerm(1.0, x[i][(n-1)][d]);
				}
				cplex.addEq(c21,1.0);
			}
			/********************/

			/********************/
			//3. Balance entre depots
			//			for(int d=0;d<inp.getVehNumber();d++) {
			//				for(int o=1; o<(n-1);o++) {
			//					IloLinearNumExpr C3 = cplex.linearNumExpr();
			//					C3.addTerm(1.0, x[0][o][d]);
			//					C3.addTerm(-1.0, x[o][n-1][d]);
			//					cplex.addEq(C3, 0); 
			//				}
			//
			//			}

			/********************/
			//4. Coneccion entre depots
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr C4 = cplex.linearNumExpr();
				C4.addTerm(1.0, x[0][(n-1)][d]);
				cplex.addEq(C4,0.0);
			}
			/********************/




			/********************/
			//6. Nada llega al depot inicial
			IloLinearNumExpr C6 = cplex.linearNumExpr();
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int i= 0; i<n;i++) {
					C6.addTerm(1.0, x[i][0][d]);
					cplex.addEq(C6, 0);
				}
			}

			/********************/

			/********************/
			//7. Nada sale del depot final
			IloLinearNumExpr C7 = cplex.linearNumExpr();
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int j= 0; j<n;j++) {
					C7.addTerm(1.0, x[n-1][j][d]);
					cplex.addEq(C7, 0);
				}
			}
			/********************/
			/********************/

			/*B. CUSTOMER*/

			/********************/ 
			//5. Un cliente puede ser visitado max una vez
			for(int j=1;j<n-1;j++) {
				IloLinearNumExpr C8 = cplex.linearNumExpr();
				for(int i=0;i<n;i++) {
					for(int d=0;d<inp.getVehNumber();d++) {
						C8.addTerm(1.0, x[i][j][d]);}}
				cplex.addLe(C8,1);
			}

			/********************/

			/********************/ 
			//5.1. Un cliente puede ser visitado max una vez
			for(int i=1;i<n-1;i++) {
				IloLinearNumExpr C81 = cplex.linearNumExpr();
				for(int j=1;j<n-1;j++) {
					for(int d=0;d<inp.getVehNumber();d++) {
						C81.addTerm(1.0, x[i][j][d]);}}
				cplex.addLe(C81,1);
			}

			/********************/


			/********************/
			//6. el numero de vehiculos que viajan en un segmento es max 1
			for(int i=1;i<n-1;i++) {
				for(int j=0;j<n-1;j++) {
					IloLinearNumExpr C9=cplex.linearNumExpr();
					for(int d=0;d<inp.getVehNumber();d++) {
						C9.addTerm(1.00,x[i][j][d]);

					}
					cplex.addLe(C9,1);
				}
			}
			/********************/


			/********************/
			//7. Para cada destino j tiene un origen i
			for(int d=0;d<inp.getVehNumber();d++){
				for(int j=1;j<n;j++){
					IloLinearNumExpr C10=cplex.linearNumExpr();
					for(int i=0;i<n-1;i++) {
						C10.addTerm(1.00,x[i][j][d]);
					}
					cplex.addLe(C10,1);}}

			/********************/


			/********************/
			//8. Para cada origen i tiene un destino j
			for(int d=0;d<inp.getVehNumber();d++){
				for(int i=0;i<n-1;i++){
					IloLinearNumExpr C11=cplex.linearNumExpr();
					for(int j=1;i<n;i++) {
						C11.addTerm(1.00,x[i][j][d]);
					}
					cplex.addLe(C11,1);}}
			/********************/



			/********************/
			//9. flujo entre rutas
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int o=1; o<(n-1);o++) {
					IloLinearNumExpr C12 = cplex.linearNumExpr();
					for(int i= 0; i<n-1;i++) {
						for(int j= 1; j<n;j++) {
							C12.addTerm(1.0, x[i][o][d]);
							C12.addTerm(-1.0, x[o][j][d]);
						}
					}
					cplex.addEq(C12, 0); 
				}

			}
			/********************/


			/********************/
			// 9. Subtour(i,j,d)$(ord(i)>1 and ord(i)<end and ord(j)>1 and ord(j)<end).. U(i,d)-U(j,d)+(card(i)*X(i,j,d))=l=card(i)-1; 
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int i=0;i<n;i++) {
					for(int j= 0; j<n;j++) {
						IloLinearNumExpr C13 = cplex.linearNumExpr();
						C13.addTerm(1.0, w[i][d]);
						C13.addTerm(-1.0, w[j][d]);
						C13.addTerm(n, x[i][j][d]);
						cplex.addLe(C13, n-1);
					}
				}
			}
			/********************/



			for(int d=0;d<inp.getVehNumber();d++) {

				IloLinearNumExpr expr10 = cplex.linearNumExpr();
				for(int i= 0; i<n;i++) {
					for(int j= 0; j<n;j++) {
						expr10.addTerm(TV[i][j], x[i][j][d]);

					}

				}
				expr10.addTerm(-1, lenght[d]);
				cplex.addEq(expr10,0);

			}

			/********************/
			//10. LenghtRoute(d)..Lenght(d)=e=sum((i,j),X(i,j,d)*(TV(i,j)));
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr C14 = cplex.linearNumExpr();
				C14.addTerm(1, lenght[d]);
				cplex.addLe(C14,inp.gettMax());

			}
			/********************/

			/*SOLVING*/
			cplex.setParam(IloCplex.Param.TimeLimit, 1000);
			if (cplex.solve()) {
				String CPLEX_file= new String(t.getInstanceName()+t.getLongSim()+"_"+"_"+t.getVariance()+"Stoch_CPLEX.txt");
				writeFile(CPLEX_file,cplex,x,V,n,K,lenght);
				//PrintWriter bw = new PrintWriter(CPLEX_file);
				System.out.println("obj = " + cplex.getObjValue());
				for (int d=0; d<inp.getVehNumber(); d++){
					System.out.println("Route "+d);
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
	static void writeFile(String filename, IloCplex cplex, IloNumVar[][][] x,int V, int n, int K, IloNumVar[] lenght) {
		try {
			PrintWriter bw = new PrintWriter(filename);
			try {
				bw.printf(Locale.US, "obj = "+ cplex.getObjValue());
				bw.println();
				for (int d=0; d<V; d++){
					bw.println("Lenght "+cplex.getValue(lenght[d]) );

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