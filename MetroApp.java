import java.util.ArrayList;
import java.util.HashMap;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.LinkedList;
import java.util.StringTokenizer;

public class MetroApp {
    public class Vertex {
        HashMap<String, Integer> nbrs = new HashMap<>();
    }
    static HashMap<String, Vertex> vtces;
    
    public MetroApp() {
        vtces = new HashMap<>();
    }
    
	public int numVetex() 
	{
		return this.vtces.size();
	}
	
	public boolean containsVertex(String vname) 
	{
		return this.vtces.containsKey(vname);
	}
	
	public void addVertex(String vname) 
	{
		Vertex vtx = new Vertex();
		vtces.put(vname, vtx);
	}

	public void removeVertex(String vname) 
	{
		Vertex vtx = vtces.get(vname);
		ArrayList<String> keys = new ArrayList<>(vtx.nbrs.keySet());

		for (String key : keys) 
		{
			Vertex nbrVtx = vtces.get(key);
			nbrVtx.nbrs.remove(vname);
		}

		vtces.remove(vname);
	}

	public int numEdges() 
	{
		ArrayList<String> keys = new ArrayList<>(vtces.keySet());
		int count = 0;

		for (String key : keys) 
		{
			Vertex vtx = vtces.get(key);
			count = count + vtx.nbrs.size();
		}

		return count / 2;
	}

	public boolean containsEdge(String vname1, String vname2) 
	{
		Vertex vtx1 = vtces.get(vname1);
		Vertex vtx2 = vtces.get(vname2);
		
		if (vtx1 == null || vtx2 == null || !vtx1.nbrs.containsKey(vname2)) {
			return false;
		}

		return true;
	}

	public void addEdge(String vname1, String vname2, int value) 
	{
		Vertex vtx1 = vtces.get(vname1); 
		Vertex vtx2 = vtces.get(vname2); 

		if (vtx1 == null || vtx2 == null || vtx1.nbrs.containsKey(vname2)) {
			return;
		}

		vtx1.nbrs.put(vname2, value);
		vtx2.nbrs.put(vname1, value);
	}

	public void removeEdge(String vname1, String vname2) 
	{
		Vertex vtx1 = vtces.get(vname1);
		Vertex vtx2 = vtces.get(vname2);
		
		//check if the vertices given or the edge between these vertices exist or not
		if (vtx1 == null || vtx2 == null || !vtx1.nbrs.containsKey(vname2)) {
			return;
		}

		vtx1.nbrs.remove(vname2);
		vtx2.nbrs.remove(vname1);
	}

    public static void CreateMetroMap(MetroApp g) {
        g.addVertex("Chennai Central");
		g.addVertex("Egmore Metro");
		g.addVertex("Nehru Park");
		g.addVertex("Kilpauk");
		g.addVertex("Pachiyappa's College");
		g.addVertex("Shenoy Nagar");
		g.addVertex("Anna Nagar East");
		g.addVertex("Anna Nagar Tower");
		g.addVertex("Thirumangalam");
		g.addVertex("Koyambedu");
		g.addVertex("CMBT");
		g.addVertex("Arumbakkam");
		g.addVertex("Vadapalani");
		g.addVertex("Ashok Nagar");
		g.addVertex("Ekkattuthangal");
		g.addVertex("Alandur");
		g.addVertex("St. Thomas Mount");
		g.addVertex("Airport");
		g.addVertex("Meenambakkam Metro");
		g.addVertex("Nanganallur Road");
		g.addVertex("Guindy");
		g.addVertex("Little Mount");
		g.addVertex("Saidapet Metro");
		g.addVertex("Nandanam");
		g.addVertex("Teynampet");
		g.addVertex("AG-DMS");
		g.addVertex("Thousand Lights");
		g.addVertex("LIC");
		g.addVertex("Government Estate");
		g.addVertex("High Court");
		g.addVertex("Mannadi");
		g.addVertex("Washermenpet");
		g.addVertex("Chennai Beach");
        g.addVertex("Chennai Fort");
        g.addVertex("Chennai Park Town");
        g.addVertex("Chintadripet");
        g.addVertex("Chepauk");
        g.addVertex("Thiruvallikeni");
        g.addVertex("Light House");
        g.addVertex("Mundagakanniamman Koil");
        g.addVertex("Thirumayilai");
        g.addVertex("Mandaveli");
        g.addVertex("Greenways Road");
        g.addVertex("Kotturpuram");
        g.addVertex("Kasturba Nagar");
        g.addVertex("Indira Nagar");
        g.addVertex("Thiruvanmiyur");
        g.addVertex("Taramani");
        g.addVertex("Perungudi");
        g.addVertex("Velachery");
        g.addVertex("Puzhuthivakkam");
        g.addVertex("Adambakkam");
        g.addVertex("Korattur");
        g.addVertex("Villivakkam");
        g.addVertex("Vyasarpadi Jeeva");
        g.addVertex("Perambur");
        g.addVertex("Perambur Carriage Works");
        g.addVertex("Perambur Loco Works");
        g.addVertex("Korukkupet");
        g.addVertex("Washermanpet");
        g.addVertex("Royapuram");
        g.addVertex("Basin Bridge");
        g.addVertex("Chennai Park");
        g.addVertex("Chetpet");
        g.addVertex("Nungambakkam");
        g.addVertex("Kodambakkam");
        g.addVertex("Mambalam");
        g.addVertex("Saidapet");
        g.addVertex("Guindy");
        g.addVertex("Pazhavanthangal");
		
		g.addEdge("Chennai Central", "Egmore Metro", 1);
		g.addEdge("Egmore Metro", "Nehru Park",1);
		g.addEdge("Nehru Park", "Kilpauk", 1);
		g.addEdge("Kilpauk", "Pachiyappa's College", 1);
		g.addEdge("Pachiyappa's College", "Shenoy Nagar", 1);
		g.addEdge("Shenoy Nagar", "Anna Nagar East", 1);
		g.addEdge("Anna Nagar East", "Anna Nagar Tower", 1);
		g.addEdge("Anna Nagar Tower", "Thirumangalam", 1);
		g.addEdge("Thirumangalam", "Koyambedu", 1);
		g.addEdge("Koyambedu", "CMBT", 1);
		g.addEdge("CMBT", "Arumbakkam",1);
		g.addEdge("Arumbakkam", "Vadapalani", 1);
		g.addEdge("Vadapalani", "Ashok Nagar", 1);
		g.addEdge("Ashok Nagar", "Ekkattuthangal", 2);
		g.addEdge("Ekkattuthangal", "Alandur", 1);
		g.addEdge("Alandur", "St. Thomas Mount", 1);
		g.addEdge("Airport", "Meenambakkam Metro", 2);
		g.addEdge("Meenambakkam Metro", "Nanganallur Road", 1);
		g.addEdge("Nanganallur Road", "Alandur", 1);
		g.addEdge("Alandur", "Guindy", 1);
		g.addEdge("Guindy", "Little Mount", 1);
		g.addEdge("Little Mount", "Saidapet Metro", 1);
		g.addEdge("Saidapet Metro", "Nandanam", 2);
		g.addEdge("Nandanam", "Teynampet", 1);
		g.addEdge("Teynampet", "AG-DMS", 1);
		g.addEdge("AG-DMS", "Thousand Lights", 1);
		g.addEdge("Thousand Lights", "LIC", 1);
		g.addEdge("LIC", "Government Estate", 1);
		g.addEdge("Government Estate", "Chennai Central", 2);
		g.addEdge("Chennai Central", "High Court", 2);
		g.addEdge("High Court", "Mannadi", 1);
		g.addEdge("Mannadi", "Washermenpet",2);
		g.addEdge("Chennai Beach", "Chennai Fort", 2);
        g.addEdge("Chennai Fort", "Chennai Park Town", 1);
        g.addEdge("Chennai Park Town", "Chintadripet", 1);
        g.addEdge("Chintadripet", "Chepauk", 2);
        g.addEdge("Chepauk", "Thiruvallikeni", 1);
        g.addEdge("Thiruvallikeni", "Light House", 1);
        g.addEdge("Light House", "Mundagakanniamman Koil", 1);
        g.addEdge("Mundagakanniamman Koil", "Thirumayilai", 1);
        g.addEdge("Thirumayilai", "Mandaveli", 1);
        g.addEdge("Mandaveli", "Greenways Road", 1);
        g.addEdge("Greenways Road", "Kotturpuram", 1);
        g.addEdge("Kotturpuram", "Kasturba Nagar", 1);
        g.addEdge("Kasturba Nagar", "Indira Nagar", 1);
        g.addEdge("Indira Nagar", "Thiruvanmiyur", 1);
        g.addEdge("Thiruvanmiyur", "Taramani", 2);
        g.addEdge("Taramani", "Perungudi", 1);
        g.addEdge("Perungudi", "Velachery", 2);
        g.addEdge("Velachery", "Puzhuthivakkam", 2);
        g.addEdge("Puzhuthivakkam", "Adambakkam", 1);
        g.addEdge("Adambakkam", "St. Thomas Mount", 1);
        g.addEdge("Korattur", "Villivakkam", 3);
        g.addEdge("Villivakkam", "Perambur Loco Works", 1);
        g.addEdge("Perambur Loco Works", "Perambur Carriage Works", 2);
        g.addEdge("Perambur Carriage Works", "Perambur", 1);
        g.addEdge("Perambur", "Vyasarpadi Jeeva", 3);
        g.addEdge("Vyasarpadi Jeeva", "Korukkupet", 3);
        g.addEdge("Vyasarpadi Jeeva", "Washermanpet", 3);
        g.addEdge("Vyasarpadi Jeeva", "Basin Bridge", 3);
        g.addEdge("Basin Bridge", "Chennai Central", 1);
        g.addEdge("Basin Bridge","Washermanpet",2);
        g.addEdge("Washermanpet", "Royapuram", 2);
        g.addEdge("Royapuram", "Chennai Beach", 2);
        g.addEdge("Chennai Beach", "Chennai Park", 3);
        g.addEdge("Chennai Park", "Chetpet", 3);
        g.addEdge("Chetpet", "Nungambakkam", 1);
        g.addEdge("Nungambakkam", "Kodambakkam", 1);
        g.addEdge("Kodambakkam", "Mambalam", 1);
        g.addEdge("Mambalam", "Saidapet", 2);
        g.addEdge("Saidapet", "Guindy", 3);
        g.addEdge("Guindy", "St. Thomas Mount", 1);
        g.addEdge("St. Thomas Mount", "Pazhavanthangal", 1);
    }
	public void display_Stations() 
	{
		System.out.println("\n***********************************************************************\n");
		ArrayList<String> keys = new ArrayList<>(vtces.keySet());
		int i=1;
		for(String key : keys) 
		{
			System.out.println(i + ". " + key);
			i++;
		}
		System.out.println("\n***********************************************************************\n");
	}
		
	public boolean hasPath(String vname1, String vname2, HashMap<String, Boolean> processed) 
	{
		
		if (containsEdge(vname1, vname2)) {
			return true;
		}

		processed.put(vname1, true);

		Vertex vtx = vtces.get(vname1);
		ArrayList<String> nbrs = new ArrayList<>(vtx.nbrs.keySet());

		for (String nbr : nbrs) 
		{

			if (!processed.containsKey(nbr))
				if (hasPath(nbr, vname2, processed))
					return true;
		}

		return false;
	}
	
	
	private class DijkstraPair implements Comparable<DijkstraPair> 
	{
		String vname;
		String psf;
		int cost;
		@Override
		public int compareTo(DijkstraPair o) 
		{
			return o.cost - this.cost;
		}
	}
	
	public int dijkstra(String src, String des, boolean nan) 
	{
		int val = 0;
		ArrayList<String> ans = new ArrayList<>();
		HashMap<String, DijkstraPair> map = new HashMap<>();

		Heap<DijkstraPair> heap = new Heap<>();

		for (String key : vtces.keySet()) 
		{
			DijkstraPair np = new DijkstraPair();
			np.vname = key;
			//np.psf = "";
			np.cost = Integer.MAX_VALUE;

			if (key.equals(src)) 
			{
				np.cost = 0;
				np.psf = key;
			}

			heap.add(np);
			map.put(key, np);
		}

		while (!heap.isEmpty()) 
		{
			DijkstraPair rp = heap.remove();
			
			if(rp.vname.equals(des))
			{
				val = rp.cost;
				break;
			}
			
			map.remove(rp.vname);

			ans.add(rp.vname);
			
			Vertex v = vtces.get(rp.vname);
			for (String nbr : v.nbrs.keySet()) 
			{
				if (map.containsKey(nbr)) 
				{
					int oc = map.get(nbr).cost;
					Vertex k = vtces.get(rp.vname);
					int nc;
					if(nan)
						nc = rp.cost + 120 + 40*k.nbrs.get(nbr);
					else
						nc = rp.cost + k.nbrs.get(nbr);

					if (nc < oc) 
					{
						DijkstraPair gp = map.get(nbr);
						gp.psf = rp.psf + nbr;
						gp.cost = nc;

						heap.updatePriority(gp);
					}
				}
			}
		}
		return val;
	}
	
	private class Pair 
	{
		String vname;
		String psf;
		int min_dis;
		int min_time;
	}
	
	public String Get_Minimum_Distance(String src, String dst) 
	{
	    int count=0;
		int min = Integer.MAX_VALUE;
		String ans = "";
		HashMap<String, Boolean> processed = new HashMap<>();
		LinkedList<Pair> stack = new LinkedList<>();
		Pair sp = new Pair();
		sp.vname = src;
		sp.psf = src + "  ";
		sp.min_dis = 0;
		sp.min_time = 0;
		stack.addFirst(sp);
		while (!stack.isEmpty()) 
		{
			Pair rp = stack.removeFirst();
			//System.out.println("VPNAME: "+rp.vname);
			//System.out.println("PSF: "+rp.vname);
			
			if (processed.containsKey(rp.vname)) 
			{
				continue;
			}
			if (!rp.vname.equals(dst)){
			processed.put(rp.vname, true);
			}
			if (rp.vname.equals(dst)) 
			{
				int temp = rp.min_dis;
				if(temp<min) {
					ans = rp.psf;
					min = temp;
				}
				if (count<10){
				processed.clear();
				count++;
				}
				continue;
			}
			Vertex rpvtx = vtces.get(rp.vname);
			ArrayList<String> nbrs = new ArrayList<>(rpvtx.nbrs.keySet());
            //System.out.println("NBRS: "+nbrs);
			for(String nbr : nbrs) 
			{
				if (!processed.containsKey(nbr)) {
					Pair np = new Pair();
					np.vname = nbr;
					np.psf = rp.psf + nbr + "  ";
					np.min_dis = rp.min_dis + rpvtx.nbrs.get(nbr); 
					if (count%2==0){
					stack.addLast(np);
					}
					else{
					    stack.addFirst(np);
					}
				}
			}
		}
		ans = ans + Integer.toString(min);
		return ans;
	}
	
	
	public String Get_Minimum_Time(String src, String dst) 
	{
		int min = Integer.MAX_VALUE;
		String ans = "";
		HashMap<String, Boolean> processed = new HashMap<>();
		LinkedList<Pair> stack = new LinkedList<>();
		Pair sp = new Pair();
		sp.vname = src;
		sp.psf = src + "  ";
		sp.min_dis = 0;
		sp.min_time = 0;
		stack.addFirst(sp);
		while (!stack.isEmpty()) {
			Pair rp = stack.removeFirst();

			if (processed.containsKey(rp.vname)) 
			{
				continue;
			}
			processed.put(rp.vname, true);
			if (rp.vname.equals(dst)) 
			{
				int temp = rp.min_time;
				if(temp<min) {
					ans = rp.psf;
					min = temp;
				}
				continue;
			}

			Vertex rpvtx = vtces.get(rp.vname);
			ArrayList<String> nbrs = new ArrayList<>(rpvtx.nbrs.keySet());

			for (String nbr : nbrs) 
			{
				if (!processed.containsKey(nbr)) {
					Pair np = new Pair();
					np.vname = nbr;
					np.psf = rp.psf + nbr + "  ";
					np.min_time = rp.min_time + 120 + 40*rpvtx.nbrs.get(nbr); 
					stack.addLast(np);
				}
			}
		}
		Double minutes = Math.ceil((double)min / 60);
		ans = ans + Double.toString(minutes);
		return ans;
	}
	
	public ArrayList<String> get_Interchanges(String str)
	{
		ArrayList<String> arr = new ArrayList<>();
		String res[] = str.split("  ");
		arr.add(res[0]);
		int count = 0;
		for(int i=1;i<res.length-1;i++)
		{
			int index = res[i].indexOf('~');
			String s = res[i].substring(index+1);
			
			if(s.length()==2)
			{
				String prev = res[i-1].substring(res[i-1].indexOf('~')+1);
				String next = res[i+1].substring(res[i+1].indexOf('~')+1);
				
				if(prev.equals(next)) 
				{
					arr.add(res[i]);
				}
				else
				{
					arr.add(res[i]+" ==> "+res[i+1]);
					i++;
					count++;
				}
			}
			else
			{
				arr.add(res[i]);
			}
		}
		arr.add(Integer.toString(count));
		arr.add(res[res.length-1]);
		return arr;
	}
	public static String[] printCodelist()
	{
		System.out.println("List of station along with their codes:\n");
		ArrayList<String> keys = new ArrayList<>(vtces.keySet());
		int i=1,j=0,m=1;
		StringTokenizer stname;
		String temp="";
		String codes[] = new String[keys.size()];
		char c;
		for(String key : keys) 
		{
			stname = new StringTokenizer(key);
			codes[i-1] = "";
			j=0;
			while (stname.hasMoreTokens())
			{
			        temp = stname.nextToken();
			        c = temp.charAt(0);
			        while (c>47 && c<58)
			        {
			                codes[i-1]+= c;
			                j++;
			                c = temp.charAt(j);
			        }
			        if ((c<48 || c>57) && c<123)
			                codes[i-1]+= c;
			}
			if (codes[i-1].length() < 2)
				codes[i-1]+= Character.toUpperCase(temp.charAt(1));
			            
			System.out.print(i + ". " + key + "\t");
			if (key.length()<(22-m))
                			System.out.print("\t");
			if (key.length()<(14-m))
                			System.out.print("\t");
                		if (key.length()<(6-m))
                			System.out.print("\t");
                		System.out.println(codes[i-1]);
			i++;
			if (i == (int)Math.pow(10,m))
			        m++;
		}
		return codes;
	}
	public class Heap<T extends Comparable<T>> {
        ArrayList<T> data = new ArrayList<>();
        HashMap<T, Integer> map = new HashMap<>();
        public void add(T item) 
    	{
    		data.add(item);   
    		map.put(item, this.data.size() - 1);
    		upheapify(data.size() - 1);
    	}

    	private void upheapify(int ci) 
    	{
    		int pi = (ci - 1) / 2;
    		if (isLarger(data.get(ci), data.get(pi)) > 0) 
    		{
    			swap(pi, ci);
    			upheapify(pi);
    		}
    	}
    
    	private void swap(int i, int j) 
    	{
    		T ith = data.get(i);
    		T jth = data.get(j);
    		
    		data.set(i, jth);
    		data.set(j, ith);
    		map.put(ith, j);
    		map.put(jth, i);
    	}
    
    	public void display() 
    	{
    		System.out.println(data);
    	}
    
    	public int size() 
    	{
    		return this.data.size();
    	}
    
    	public boolean isEmpty() 
    	{
    		return this.size() == 0;
    	}
    
    	public T remove() 
    	{
    		swap(0, this.data.size() - 1);
    		T rv = this.data.remove(this.data.size() - 1);
    		downheapify(0);
    
    		map.remove(rv);
    		return rv;
    	}
    
    	private void downheapify(int pi) 
    	{
    		int lci = 2 * pi + 1;
    		int rci = 2 * pi + 2;
    		int mini = pi;
    
    		if (lci < this.data.size() && isLarger(data.get(lci), data.get(mini)) > 0)
    		{
    			mini = lci;
    		}
    		
    		if (rci < this.data.size() && isLarger(data.get(rci), data.get(mini)) > 0) 
    		{
    			mini = rci;
    		}
    		
    		if (mini != pi)
    		{
    			swap(mini, pi);
    			downheapify(mini);
    		}
    	}
    
    	public T get() 
    	{
    		return this.data.get(0);
    	}
    
    	public int isLarger(T t, T o) 
    	{
    		return t.compareTo(o);
    	}
    
    	public void updatePriority(T pair) 
    	{
    		int index = map.get(pair);
    		upheapify(index);
    	}
    }
    void displayMap()
    {
        System.out.println("Line-1 (Green Map) ");
        System.out.println("--------------------");
        System.out.println("Central Metro");
        System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Egmore Metro");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Nehru Park");
		System.out.println("       |");
        System.out.println("       ↓");
	    System.out.println("Kilpauk");
	    System.out.println("       |");
        System.out.println("       ↓");
        System.out.println("Pachiyappa's College");
        System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Shenoy Nagar");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Anna Nagar East");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Anna Nagar Tower");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Thirumangalam");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Koyambedu");
		System.out.println("       |");
        System.out.println("       ↓");
        System.out.println("CMBT");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Arumbakkam");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Vadapalani");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Ashok Nagar");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Ekkattuthangal");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Alandur");
		System.out.println("       |");
        System.out.println("       ↓");
	    System.out.println("St. Thomas Mount Metro");
	    System.out.println("");
	    System.out.println("************************************");
	    System.out.println("");
	    System.out.println("Line-2 (Red Line)");
	    System.out.println("-----------------");
	    System.out.println("Airport");
	    System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Meenambakkam Metro");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Nanganallur Road");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Guindy");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Little Mount");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Saidapet Metro");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Nandanam");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Teynampet");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("AG-DMS");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Thousand Lights");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("LIC");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Government Estate");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("High Court");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Mannadi");
		System.out.println("       |");
        System.out.println("       ↓");
		System.out.println("Washermenpet Metro");
		System.out.println("");
		System.out.println("****************************");
		System.out.println("");
		System.out.println("Line-3 (MRTS)");
		System.out.println("--------------");
        System.out.println("Chennai Beach");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Chennai Fort");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Chennai Park Town");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Chintadripet");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Chepauk");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Thiruvallikeni");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Light House");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Mundagakanniamman Koil");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Thirumayilai");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Mandaveli");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Greenways Road");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Kotturpuram");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Kasturba Nagar");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Indira Nagar");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Thiruvanmiyur");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Taramani");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Perungudi");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Velachery");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Puzhuthivakkam");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Adambakkam");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("St Thomas Mount");
        System.out.println("");
        System.out.println("***********************");
        System.out.println("");
        System.out.println("Line-4 (SOUTHERN RAILWAY)");
        System.out.println("-------------------------");
        System.out.println("Korattur");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Villivakkam");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Vyasarpadi Jeeva");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Perambur");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Perambur Carriage Works");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Perambur Loco Works");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Korukkupet");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Washermanpet");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Royapuram");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Basin Bridge");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Chennai Park");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Chetpet");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Nungambakkam");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Kodambakkam");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Mambalam");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Saidapet");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Guindy");
        System.out.println("   |");
        System.out.println("   ↓");
        System.out.println("Pazhavanthangal");
        System.out.println("****************************");
    }
       
	
    public static void main(String[] args) throws IOException {
        MetroApp g = new MetroApp();
        CreateMetroMap(g);

        System.out.println("\n\t\t\t****WELCOME TO THE METRO APP*****");
        BufferedReader inp = new BufferedReader(new InputStreamReader(System.in));
        String[] metrogreen = {"Central Metro","Egmore Metro","Nehru Park","Kilpauk","Pachiyappa's College","Shenoy Nagar","Anna Nagar East","Anna Nagar Tower","Thirumangalam","Koyambedu","CMBT","Arumbakkam","Vadapalani","Ashok Nagar","Ekkattuthangal","Alandur","St. Thomas Mount Metro"};
        String[] mrts = {
            "Chennai Beach", "Chennai Fort", "Chennai Park Town", "Chintadripet",
            "Chepauk", "Thiruvallikeni", "Light House", "Mundagakanniamman Koil",
            "Thirumayilai", "Mandaveli", "Greenways Road", "Kotturpuram",
            "Kasturba Nagar", "Indira Nagar", "Thiruvanmiyur", "Taramani",
            "Perungudi", "Velachery", "Puzhuthivakkam", "Adambakkam",
            "St. Thomas Mount"
        };
        String[] southern = {
            "Korattur", "Villivakkam", "Perambur Loco Works", "Perambur Carriage Works",
            "Perambur", "Vyasarpadi Jeeva", "Korukkupet", "Washermanpet", "Basin Bridge",
            "Chennai Central", "Royapuram", "Chennai Beach", "Chennai Park", "Chetpet",
            "Nungambakkam", "Kodambakkam", "Mambalam", "Saidapet", "Guindy", "St. Thomas Mount",
            "Pazhavanthangal"
        };
        String[] metrored = {"Airport","Meenambakkam Metro","Nanganallur Road","Alandur","Guindy","Little Mount","Saidapet Metro","Nandanam","Teynampet","AG-DMS","Thousand Lights","LIC","Government Estate","High Court","Mannadi","Washermenpet Metro","Central Metro"};
        while (true) {
            System.out.println("\t\t\t\t~~LIST OF ACTIONS~~\n\n");
			System.out.println("1. LIST ALL THE STATIONS IN THE MAP");
			System.out.println("2. SHOW THE METRO MAP");
			System.out.println("3. GET SHORTEST DISTANCE FROM A 'SOURCE' STATION TO 'DESTINATION' STATION");
			System.out.println("4. GET SHORTEST TIME TO REACH FROM A 'SOURCE' STATION TO 'DESTINATION' STATION");
			System.out.println("5. GET SHORTEST PATH (DISTANCE WISE) TO REACH FROM A 'SOURCE' STATION TO 'DESTINATION' STATION");
			System.out.println("6. GET SHORTEST PATH (TIME WISE) TO REACH FROM A 'SOURCE' STATION TO 'DESTINATION' STATION");
			System.out.println("7. EXIT THE MENU");
			System.out.print("\nENTER YOUR CHOICE FROM THE ABOVE LIST (1 to 7) : ");
			int choice = -1;
			try {
				choice = Integer.parseInt(inp.readLine());
			} catch(Exception e) {
	
			}
			System.out.print("\n***********************************************************\n");
			if(choice == 7)
			{
				System.exit(0);
			}
			switch(choice)
			{
			case 1:
				g.display_Stations();
				break;
		
			case 2:
				g.displayMap();
				break;
			
			case 3:
				ArrayList<String> keys = new ArrayList<>(vtces.keySet());
				String codes[] = printCodelist();
				System.out.println("\n1. TO ENTER SERIAL NO. OF STATIONS\n2. TO ENTER CODE OF STATIONS\n3. TO ENTER NAME OF STATIONS\n");
				System.out.println("ENTER YOUR CHOICE:");
			        int ch = Integer.parseInt(inp.readLine());
				int j;
					
				String st1 = "", st2 = "";
				System.out.println("ENTER THE SOURCE AND DESTINATION STATIONS");
				if (ch == 1)
				{
				    st1 = keys.get(Integer.parseInt(inp.readLine())-1);
				    st2 = keys.get(Integer.parseInt(inp.readLine())-1);
				}
				else if (ch == 2)
				{
				    String a,b;
				    a = (inp.readLine()).toUpperCase();
				    for (j=0;j<keys.size();j++)
				       if (a.equals(codes[j]))
				           break;
				    st1 = keys.get(j);
				    b = (inp.readLine()).toUpperCase();
				    for (j=0;j<keys.size();j++)
				       if (b.equals(codes[j]))
				           break;
				    st2 = keys.get(j);
				}
				else if (ch == 3)
				{
				    st1 = inp.readLine();
				    st2 = inp.readLine();
				}
				else
				{
				    System.out.println("Invalid choice");
				    System.exit(0);
				}
			
				HashMap<String, Boolean> processed = new HashMap<>();
				if(!g.containsVertex(st1) || !g.containsVertex(st2) || !g.hasPath(st1, st2, processed))
					System.out.println("THE INPUTS ARE INVALID");
				else
				System.out.println("SHORTEST DISTANCE FROM "+st1+" TO "+st2+" IS "+g.dijkstra(st1, st2, false)+"KM\n");
				break;
			
			case 4:
                ArrayList<String> key = new ArrayList<>(vtces.keySet());
				String code[] = printCodelist();
				System.out.println("\n1. TO ENTER SERIAL NO. OF STATIONS\n2. TO ENTER CODE OF STATIONS\n3. TO ENTER NAME OF STATIONS\n");
				System.out.println("ENTER YOUR CHOICE:");
			    int cha = Integer.parseInt(inp.readLine());
				int ja;
				String sat1 = "", sat2 = "";
				System.out.println("ENTER THE SOURCE AND DESTINATION STATIONS");
				if (cha == 1)
				{
				    sat1 = key.get(Integer.parseInt(inp.readLine())-1);
				    sat2 = key.get(Integer.parseInt(inp.readLine())-1);
				}
				else if (cha == 2)
				{
				    String a,b;
				    a = (inp.readLine()).toUpperCase();
				    for (ja=0;ja<key.size();ja++)
				       if (a.equals(code[ja]))
				           break;
				    sat1 = key.get(ja);
				    b = (inp.readLine()).toUpperCase();
				    for (ja=0;ja<key.size();ja++)
				       if (b.equals(code[ja]))
				           break;
				    sat2 = key.get(ja);
				}
				else if (cha == 3)
				{
				    sat1 = inp.readLine();
				    sat2 = inp.readLine();
				}
				else
				{
				    System.out.println("Invalid choice");
				    System.exit(0);
				}
			
				HashMap<String, Boolean> processed1= new HashMap<>();				
				System.out.println("SHORTEST TIME FROM ("+sat1+") TO ("+sat2+") IS "+g.dijkstra(sat1, sat2, true)/60+" MINUTES\n\n");
				break;
			
			case 5:
                ArrayList<String> keys1 = new ArrayList<>(vtces.keySet());
				String codes1[] = printCodelist();
				System.out.println("\n1. TO ENTER SERIAL NO. OF STATIONS\n2. TO ENTER CODE OF STATIONS\n3. TO ENTER NAME OF STATIONS\n");
				System.out.println("ENTER YOUR CHOICE:");
			    int ch1 = Integer.parseInt(inp.readLine());
				int j1;
				String s1 = "", s2 = "";
				System.out.println("ENTER THE SOURCE AND DESTINATION STATIONS");
				if (ch1 == 1)
				{
				    s1 = keys1.get(Integer.parseInt(inp.readLine())-1);
				    s2 = keys1.get(Integer.parseInt(inp.readLine())-1);
				}
				else if (ch1 == 2)
				{
				    String a,b;
				    a = (inp.readLine()).toUpperCase();
				    for (j1=0;j1<keys1.size();j1++)
				       if (a.equals(codes1[j1]))
				           break;
				    s1 = keys1.get(j1);
				    b = (inp.readLine()).toUpperCase();
				    for (j1=0;j1<keys1.size();j1++)
				       if (b.equals(codes1[j1]))
				           break;
				    s2 = keys1.get(j1);
				}
				else if (ch1 == 3)
				{
				    s1 = inp.readLine();
				    s2 = inp.readLine();
				}
				else
				{
				    System.out.println("Invalid choice");
				    System.exit(0);
				}
			
				HashMap<String, Boolean> processed2 = new HashMap<>();
				if(!g.containsVertex(s1) || !g.containsVertex(s2) || !g.hasPath(s1, s2, processed2))
					System.out.println("THE INPUTS ARE INVALID");
				else 
				{
					ArrayList<String> str = g.get_Interchanges(g.Get_Minimum_Distance(s1, s2));
					//System.out.println(g.Get_Minimum_Distance(s1, s2));
					//System.out.println(str);
					int len = str.size();
					System.out.println("SOURCE STATION : " + s1);
					System.out.println("DESTINATION STATION : " + s2);
					System.out.println("DISTANCE : " + str.get(len-1) + " KM");
					//System.out.println("NUMBER OF INTERCHANGES : " + str.get(len-2));
					System.out.println("-----------------START-----------------");
					System.out.println(str.get(0));
					int inter,gr=0,re=0,uiy=0,qwe=0,yiu=0,ewq=0,bnm=0,mnb=0,zxc=0,cxz=0;
					int grw=0,rew=0,mew=0,poi=0;;
					for (String station : metrogreen){
					    if (station.equals(str.get(0))){
					        grw=1;
					    }
					}
					for (String station : metrored){
					    if (station.equals(str.get(0))){
					        rew=1;
					    }
					}
					for (String station : mrts) {
                        if (station.equals(str.get(0))) {
                            mew=1;
                        }
                    }
                    for (String station : southern){
					    if (station.equals(str.get(0))){
					        poi=1;
					    }
					}
					for(int i=1; i<len-2; i++)
					{
				        yiu=uiy;
					    ewq=qwe;
					    mnb=bnm;
					    cxz=zxc;
					    uiy=0;
					    qwe=0;
					    bnm=0;
					    zxc=0;
						for (String station : metrogreen) {
                            if (station.equals(str.get(i))) {
                                uiy=1;
                            }
                        }
                        for (String station1 : metrored) {
                            if (station1.equals(str.get(i))) {
                                qwe=1;
                            }
                        }
                        for (String station2 : mrts) {
                            if (station2.equals(str.get(i))) {
                                bnm=1;
                            }
                        }
                        for (String station3 : southern) {
                            if (station3.equals(str.get(i))) {
                                zxc=1;
                            }
                        }
                        
                        if (yiu==1 && ewq==1 && uiy==1 && grw==0 || yiu==1 && mnb==1 && uiy==1 &&grw==00 || yiu==1 && cxz ==1 && uiy==1 && grw==0){
                            re+=1;
                            System.out.println("---------------------------");
                            System.out.println("PLEASE CHANGE TO GREEN LINE");
                            System.out.println("---------------------------");
                            grw=1;
                        }
                        if (yiu==1 && ewq==1 && qwe==1 && rew==0 || ewq==1 && mnb==1 && qwe==1 && rew==0||ewq==1 && cxz==1 && rew==0 && qwe==1){
                            re+=1;
                            System.out.println("-------------------------");
                            System.out.println("PLEASE CHANGE TO RED LINE");
                            System.out.println("-------------------------");
                            rew=1;
                        }
                        if (bnm==1 && yiu==1 && mnb==1 && mew==0 || bnm==1 && mnb==1 && ewq==1 && mew==0||bnm==1 && mnb==1 && cxz==1 && mew==0){
                            re+=1;
                            System.out.println("-------------------------");
                            System.out.println("PLEASE CHANGE TO MRTS LINE");
                            System.out.println("-------------------------");
                            mew=1;
                        }
                        if (zxc==1 && cxz==1 && yiu==1 && poi==0 || zxc==1 && cxz==1 && ewq==1 && poi==0||cxz==1 && zxc==1 && poi==0 && mnb==0){
                            re+=1;
                            System.out.println("--------------------------------------");
                            System.out.println("PLEASE CHANGE TO SOUTHERN RAILWAY LINE");
                            System.out.println("--------------------------------------");
                            poi=1;
                        }
                        
                        System.out.println(str.get(i));
					}
					int tot=gr+re-1;
					if (tot<0){
					    tot=0;
					}
					System.out.print("----------------END---------------------");
					System.out.println("\nNUMBER OF INTERCHANGES : " + re);
					System.out.println("\n~~~~~~~~~~~~~");
				}
				break;
			
			case 6:
                ArrayList<String> keys2 = new ArrayList<>(vtces.keySet());
				String codes2[] = printCodelist();
				System.out.println("\n1. TO ENTER SERIAL NO. OF STATIONS\n2. TO ENTER CODE OF STATIONS\n3. TO ENTER NAME OF STATIONS\n");
				System.out.println("ENTER YOUR CHOICE:");
			    int ch2 = Integer.parseInt(inp.readLine());
				int j2;
				String ss1 = "", ss2 = "";
				System.out.println("ENTER THE SOURCE AND DESTINATION STATIONS");
				if (ch2 == 1)
				{
				    ss1 = keys2.get(Integer.parseInt(inp.readLine())-1);
				    ss2 = keys2.get(Integer.parseInt(inp.readLine())-1);
				}
				else if (ch2 == 2)
				{
				    String a,b;
				    a = (inp.readLine()).toUpperCase();
				    for (j2=0;j2<keys2.size();j2++)
				       if (a.equals(codes2[j2]))
				           break;
				    ss1 = keys2.get(j2);
				    b = (inp.readLine()).toUpperCase();
				    for (j2=0;j2<keys2.size();j2++)
				       if (b.equals(codes2[j2]))
				           break;
				    ss2 = keys2.get(j2);
				}
				else if (ch2 == 3)
				{
				    ss1 = inp.readLine();
				    ss2 = inp.readLine();
				}
				else
				{
				    System.out.println("Invalid choice");
				    System.exit(0);
				}
			
				HashMap<String, Boolean> processed3 = new HashMap<>();
				if(!g.containsVertex(ss1) || !g.containsVertex(ss2) || !g.hasPath(ss1, ss2, processed3))
					System.out.println("THE INPUTS ARE INVALID");
				else
				{
					ArrayList<String> str = g.get_Interchanges(g.Get_Minimum_Time(ss1, ss2));
					int len = str.size();
					System.out.println("SOURCE STATION : " + ss1);
					System.out.println("DESTINATION STATION : " + ss2);
					System.out.println("TIME : " + str.get(len-1)+" MINUTES");
					System.out.println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
                    System.out.println("-----------------START-----------------");
					System.out.println(str.get(0));
					int inter,gr=0,re=0,uiy=0,qwe=0,yiu=0,ewq=0,bnm=0,mnb=0,zxc=0,cxz=0;
					int grw=0,rew=0,mew=0,poi=0;;
					for (String station : metrogreen){
					    if (station.equals(str.get(0))){
					        grw=1;
					    }
					}
					for (String station : metrored){
					    if (station.equals(str.get(0))){
					        rew=1;
					    }
					}
					for (String station : mrts) {
                        if (station.equals(str.get(0))) {
                            mew=1;
                        }
                    }
                    for (String station : southern){
					    if (station.equals(str.get(0))){
					        poi=1;
					    }
					}
					for(int i=1; i<len-2; i++)
					{
                        yiu=uiy;
					    ewq=qwe;
					    mnb=bnm;
					    cxz=zxc;
					    uiy=0;
					    qwe=0;
					    bnm=0;
					    zxc=0;
						for (String station : metrogreen) {
                            if (station.equals(str.get(i))) {
                                uiy=1;
                            }
                        }
                        for (String station1 : metrored) {
                            if (station1.equals(str.get(i))) {
                                qwe=1;
                            }
                        }
                        for (String station2 : mrts) {
                            if (station2.equals(str.get(i))) {
                                bnm=1;
                            }
                        }
                        for (String station3 : southern) {
                            if (station3.equals(str.get(i))) {
                                zxc=1;
                            }
                        }
                        
                        if (yiu==1 && ewq==1 && uiy==1 && grw==0 || yiu==1 && mnb==1 && uiy==1 &&grw==00 || yiu==1 && cxz ==1 && uiy==1 && grw==0){
                            re+=1;
                            System.out.println("---------------------------");
                            System.out.println("PLEASE CHANGE TO GREEN LINE");
                            System.out.println("---------------------------");
                            grw=1;
                        }
                        if (yiu==1 && ewq==1 && qwe==1 && rew==0 || ewq==1 && mnb==1 && qwe==1 && rew==0||ewq==1 && cxz==1 && rew==0 && qwe==1){
                            re+=1;
                            System.out.println("-------------------------");
                            System.out.println("PLEASE CHANGE TO RED LINE");
                            System.out.println("-------------------------");
                            rew=1;
                        }
                        if (bnm==1 && yiu==1 && mnb==1 && mew==0 || bnm==1 && mnb==1 && ewq==1 && mew==0||bnm==1 && mnb==1 && cxz==1 && mew==0){
                            re+=1;
                            System.out.println("-------------------------");
                            System.out.println("PLEASE CHANGE TO MRTS LINE");
                            System.out.println("-------------------------");
                            mew=1;
                        }
                        if (zxc==1 && cxz==1 && yiu==1 && poi==0 || zxc==1 && cxz==1 && ewq==1 && poi==0||cxz==1 && zxc==1 && poi==0 && mnb==0){
                            re+=1;
                            System.out.println("--------------------------------------");
                            System.out.println("PLEASE CHANGE TO SOUTHERN RAILWAY LINE");
                            System.out.println("--------------------------------------");
                            poi=1;
                        }
                        
                        System.out.println(str.get(i));
					}
					System.out.print(str.get(len-3) + "   ==>    END");
					int tot=gr+re-1;
					if (tot<0){
					    tot=0;
					}
					System.out.println("\nNUMBER OF INTERCHANGES : " + re);
					System.out.println("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~");
				}
				break;	
           	         default: 
                	        System.out.println("Please enter a valid option! ");
                    	    System.out.println("The options you can choose are from 1 to 6. ");
                        
			}
		}
    }
}
        
