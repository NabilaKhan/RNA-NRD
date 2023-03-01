
public class Atom {
	int sn;
	Residue res;
	String atom;
	Point coord;
	
	public Atom(int sn, Residue res, String atom, Point coord){
		this.sn=sn;
		this.res=res;
		this.atom=atom;
		this.coord=coord;
	}
	
	public String toString(){
		return sn+"\t"+res.toString()+"\t"+atom+':'+coord;
	}

}
