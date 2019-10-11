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

public class TOP_StochasticVersion {

	public static void solveMe(Inputs inp, Test t) {
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/*SETS*/

		int V=inp.getVehNumber(); // set of vehicles
		int K=t.getLongSim()+1; // set of scenarios/realizations
		int n=inp.getNodes().length; //set of node
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
			U[i]=IDnode.getProfit();
		}

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

			//Ideal reward per route in each scenario
			IloNumVar[][] R=new IloNumVar[V][K];
			for(int d=0;d<V;d++){
				R[d]=cplex.numVarArray(K, 0,Double.MAX_VALUE);
			}

			//Real reward per route in each scenario
			IloNumVar[][] Reward=new IloNumVar[V][K];
			for(int d=0;d<V;d++){
				Reward[d]=cplex.numVarArray(K, 0,Double.MAX_VALUE);
			}

			/*OBJECTIVE FUNCTION*/

			IloLinearNumExpr obj = cplex.linearNumExpr();
			for(int d=0; d<V;d++) {
				for(int k=0;k<K;k++ ) {
					obj.addTerm(1, Reward[d][k]);
				}
			}
			cplex.addMaximize(obj);

			/*CONSTRAINTS*/

			//1. cust_depot(d)..sum((i,j)$(ord(j)=end),X(i,j,d))=l=1;
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr expr1 = cplex.linearNumExpr();
				for(int i= 0; i<n;i++) {
					expr1.addTerm(1.0, x[i][(n-1)][d]);
				}
				cplex.addEq(expr1,1.0);
			}

			//2.origen_cust(d)..sum((j),X('0',j,d))=l=1;  
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr expr2 = cplex.linearNumExpr();
				for(int j= 0; j<n;j++) {
					expr2.addTerm(1.0, x[0][j][d]);
				}
				cplex.addEq(expr2,1.0);
			}

			//3. origen_end(j,d)$(ord(j)=end)..,X('0',j,d)=e=0; 
			for(int d=0;d<inp.getVehNumber();d++) {
				IloLinearNumExpr expr3 = cplex.linearNumExpr();
				expr3.addTerm(1.0, x[0][(n-1)][d]);
				cplex.addEq(expr3,0.0);
			}

			//4. visit_cust(j)$(ord(j)<end and ord(j)>0)..sum((i,d),X(i,j,d))=l=1; 
			for(int j=1;j<n-1;j++) {
				IloLinearNumExpr expr4 = cplex.linearNumExpr();
				for(int i=0;i<n;i++) {
					for(int d=0;d<V;d++) {
						expr4.addTerm(1.0, x[i][j][d]);}}
				cplex.addLe(expr4,1);
			}

			//5. flow(o,d)$(ord(o)<>0 and ord(o)<>end)..sum(i,X(i,o,d))=e=sum(j,X(o,j,d)); 
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int o=1; o<(n-1);o++) {
					IloLinearNumExpr expr5 = cplex.linearNumExpr();
					for(int i= 0; i<n-1;i++) {
						for(int j= 1; j<n;j++) {
							expr5.addTerm(1.0, x[i][o][d]);
							expr5.addTerm(-1.0, x[o][j][d]);
						}
					}
					cplex.addEq(expr5, 0); 
				}

			}

			//6. star_route(i,j,d)$(ord(j)=0)..X(i,j,d)=e=0;                  
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int i= 0; i<n;i++) {
					IloLinearNumExpr expr6 = cplex.linearNumExpr();
					expr6.addTerm(1.0, x[i][0][d]);
					cplex.addEq(expr6, 0.0);
				}
			}

			//7. segment_origen(d,i)$(ord(i)<>end)..sum(j$(ord(j)<>0),X(i,j,d))=l=1;        
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int i= 0; i<(n-1);i++) {
					IloLinearNumExpr expr7 = cplex.linearNumExpr();
					for(int j= 1; j<n;j++) {
						expr7.addTerm(1.0, x[i][j][d]);
					}
					cplex.addLe(expr7, 1);
				}
			}

			//8. segment_end(d,j)$(ord(j)<>0)..sum(i$(ord(i)<>end),X(i,j,d))=l=1; 
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int j= 1; j<n;j++) {
					IloLinearNumExpr expr8 = cplex.linearNumExpr();
					for(int i= 0; i<(n-1);i++) {
						expr8.addTerm(1.0, x[i][j][d]);
					}
					cplex.addLe(expr8, 1);
				}
			}

			// 9. Subtour(i,j,d).. U(i,d)-U(j,d)+(card(i)*X(i,j,d))=l=card(i)-1; 
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int i=0;i<n;i++) {
					for(int j= 0; j<n;j++) {
						IloLinearNumExpr expr9 = cplex.linearNumExpr();
						expr9.addTerm(1.0, w[i][d]);
						expr9.addTerm(-1.0, w[j][d]);
						expr9.addTerm(n, x[i][j][d]);
						cplex.addLe(expr9, n-1);
					}
				}
			}

			//10. LenghtRoute(d,k)..Lenght(d,k)=e=sum((i,j),X(i,j,d)*(TV(i,j)));
			for(int d=0;d<inp.getVehNumber();d++) {
				for(int k=0;k<K;k++) {
					IloLinearNumExpr expr10 = cplex.linearNumExpr();
					for(int i= 0; i<n;i++) {
						for(int j= 0; j<n;j++) {
							expr10.addTerm(TV[i][j][k], x[i][j][d]);
						}
					}
					expr10.addTerm(-1, lenght[d][k]);
					cplex.addLe(expr10,0);
				}
			}

			//11. rewardperRoute(d,k)..R(d,k)=e=sum((i,j),X(i,j,d)*(U(i)));
			for(int k=0;k<K;k++) {
				for(int d=0;d<inp.getVehNumber();d++) {
					IloLinearNumExpr expr11 = cplex.linearNumExpr();
					for(int i= 0; i<n;i++) {
						for(int j= 0; j<n;j++) {
							expr11.addTerm(U[i], x[i][j][d]);
						}
					}
					expr11.addTerm(-1, R[d][k]);
					cplex.addGe(expr11,0);
				}
			}

			//12*. rewardperRoute(k)$(TV[i][j][k]*x[i][j][d]<=Tlim)..Reward(d,k)=e=sum((i,j),X(i,j,d)*(U(i)));
			for(int k=0;k<K;k++) {
				for(int d=0;d<inp.getVehNumber();d++) {
					cplex.add(cplex.ifThen(cplex.ge(lenght[d][k], inp.gettMax()), cplex.eq(Reward[d][k], 0.0)));
					cplex.add(cplex.ifThen(cplex.le(lenght[d][k], inp.gettMax()), cplex.eq(Reward[d][k],R[d][k])));
				}
			}

			/*SOLVING*/

			//STOP CRITERIA                         
			cplex.setParam(IloCplex.Param.TimeLimit, 5000);

			//SOLVER
			if (cplex.solve()) {
				String CPLEX_file= new String(t.getInstanceName()+t.getLongSim()+"_"+"_"+t.getVariance()+"Stoch_CPLEX.txt");
				writeFile(CPLEX_file,cplex,x,Reward,lenght,V,n,K);
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
			}
			cplex.end(); 
		}

		catch(IloException e) {
			e.printStackTrace();
		}

		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	}

	static void writeFile(String filename, IloCplex cplex, IloNumVar[][][] x,IloNumVar[][] R, IloNumVar[][] lenght,int V, int n, int K) {
		//CPLEX_file,cplex,x,R,lenght,V,n,K
		try {
			PrintWriter bw = new PrintWriter(filename);
			try {
				bw.printf(Locale.US, "obj = "+ cplex.getObjValue()/K);
				for (int d=0; d<V; d++){
					for (int k=0; k<K; k++){
						if(cplex.getValue(lenght[d][k])!=0) {
							bw.println();
							bw.printf(Locale.US, "lenght = "+ cplex.getValue(lenght[d][k]));}
						if(cplex.getValue(R[d][k])!=0) {
							bw.println();
							bw.printf(Locale.US, "Reward = "+ cplex.getValue(R[d][k]));}}}
				bw.println();
				for (int d=0; d<V; d++){
					bw.println();
					bw.println("Route "+d+1);
					for (int i=0; i<n; i++){
						for (int j=0; j<n; j++){
							if (i!=j && cplex.getValue(x[i][j][d])>0.9999) {
								bw.println(i+" - "+j);
							}
						}
					}
				}
				bw.println();
				bw.print(cplex.getCplexTime());
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