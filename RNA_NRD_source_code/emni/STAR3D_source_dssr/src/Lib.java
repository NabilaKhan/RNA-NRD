import java.io.*;
import java.util.*;

public class Lib {
	//get the base pair (index based) on a sequence coming from a PDB file
	public static HashSet<Pair<Integer, String>> get_seq_bp(PDBParser pdb, HashSet<Pair<ResID, String>> basepair, String chainID){
		ArrayList<ResID> rid_list=new ArrayList<ResID>();
		for(Residue R: pdb.get_chain_res().get(chainID)) {
			rid_list.add(R.rid);
		}

		HashSet<Pair<Integer, String>> bp=new HashSet<Pair<Integer, String>>();

		int i=0, j=0;
		for(Pair<ResID, String> P: basepair){
			if(P.left.chainID.equals(chainID) && P.right.chainID.equals(chainID)){
				i=rid_list.indexOf(P.left);
				j=rid_list.indexOf(P.right);
				if(i != -1 && j != -1)
					bp.add(new Pair<Integer, String>(i, j, P.v));
			}
		}
		return bp;
	}

	//get the pairing partners of bases on a sequence coming from a PDB file
	public static HashMap<Integer, ArrayList<Integer>> get_seq_pairing_map(HashSet<Pair<Integer, String>> bp){
		HashMap<Integer, ArrayList<Integer>> pairing_map=new HashMap<Integer, ArrayList<Integer>>();
		for(Pair<Integer, String> P: bp){
			if(! "WHS".contains(P.v.substring(0,1)) || ! "WHS".contains(P.v.substring(1,2))) continue;
			if(pairing_map.get(P.left)==null) pairing_map.put(P.left, new ArrayList<Integer>());
			if(pairing_map.get(P.right)==null) pairing_map.put(P.right, new ArrayList<Integer>());
			pairing_map.get(P.left).add(P.right);
			pairing_map.get(P.right).add(P.left);
		}
		return pairing_map;
	}	
	
	public static List<Pair<Integer, Integer>> get_ct_bp(File ct_fn) throws IOException{
		List<Pair<Integer, Integer>> BP=new ArrayList<Pair<Integer, Integer>>();
		BufferedReader in=new BufferedReader(new FileReader(ct_fn));
		String line;
		String[] token;
		
		in.readLine();
		int i, j;
		while((line=in.readLine())!=null){
			token=line.trim().split("\\s+");
			i=Integer.valueOf(token[0]);
			j=Integer.valueOf(token[4]);
			if(i<j) BP.add(new Pair<Integer, Integer>(i-1, j-1, 1));
		}
		in.close();
		return BP;
	}
	
	//find the neighbor of i in graph G
	private static Set<Integer> GN(int[][] G, int i){
		Set<Integer> ret=new HashSet<Integer>();
		for(int j=0; j<G[i].length; j++){
			if(G[i][j]==1 && j!=i) ret.add(j);
		}
		return ret;
	}
	
	public static void Bron_Kerbosch(int[][] G, Set<Integer> R, Set<Integer> P, Set<Integer> X, List<Set<Integer>> ret){
		if(P.isEmpty() && X.isEmpty()){
			ret.add(R);
			return;
		}
		Set<Integer> PX = new HashSet<Integer>(P);
		PX.addAll(X);
		
		int i=new Random().nextInt(PX.size());
		int j=0, u=0;
		for(Integer I: PX){
			if(j==i) u=I;
			j++;
		}
		
		Set<Integer> PGNU=new HashSet<Integer>(P);

		PGNU.removeAll(GN(G, u));
		for(Integer v: PGNU){
			Set<Integer> RV=new HashSet<Integer>(R);
			RV.add((Integer)v);
			
			Set<Integer> PGNV=new HashSet<Integer>();
			for(Integer I: GN(G, v))
				if(P.contains(I)) PGNV.add(I);
			
			Set<Integer> XGNV=new HashSet<Integer>();
			for(Integer I: GN(G, v))
				if(X.contains(I)) XGNV.add(I);
			
			Bron_Kerbosch(G, RV, PGNV, XGNV, ret);
			P.remove((Integer)v);
			X.add((Integer)v);
		}
	}
}
