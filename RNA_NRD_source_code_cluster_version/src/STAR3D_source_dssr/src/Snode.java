import java.util.*;

public class Snode {
	Pair<Integer, Integer> stack;
	List<Snode> children;
	
	Snode(Pair<Integer, Integer> stack, List<Snode> children){
		this.stack=stack; this.children=children;
	}

	public void show(int lvl){
		String head="";
		for(int i=0; i<lvl; i++) head+="+";
		System.out.println(head+stack);
		for(Snode SN: children) SN.show(lvl+1);
	}
	
	public static Snode stack_tree(List<Pair<Integer, Integer>> stack_list, int len){
		LinkedList<Snode> node_stack=new LinkedList<Snode>();
		
		Snode root=new Snode(new Pair<Integer, Integer>(-1, len, 1), new ArrayList<Snode>());
		node_stack.addFirst(root);
		for(Pair<Integer, Integer> S: stack_list){
			Snode N=new Snode(S, new ArrayList<Snode>());
			Snode M=node_stack.getFirst(); 
			while(! (M.stack.left < N.stack.left && N.stack.right < M.stack.right)){
				node_stack.removeFirst();
				M=node_stack.getFirst();
			}
			node_stack.get(0).children.add(N);
			node_stack.addFirst(N);
		}
		return node_stack.get(node_stack.size()-1);
	}
	
	public static void order_loop(Snode root, List<Pair<Integer, Integer>> ret){
		int i=root.stack.left+root.stack.v;
		int j=root.stack.right-root.stack.v;
		
		int l, k;
		for(Snode SN: root.children){
			l=SN.stack.left;
			k=SN.stack.right;
			if(i<=l-1) ret.add(new Pair<Integer, Integer>(i, l-1, 0));
			else ret.add(new Pair<Integer, Integer>(-1, -1, 0));
			i=k+1;
		}
		if(i<=j) ret.add(new Pair<Integer, Integer>(i, j, 0));
		else ret.add(new Pair<Integer, Integer>(-1, -1, 0));
		
		for(Snode SN: root.children) order_loop(SN, ret);
	}
}
