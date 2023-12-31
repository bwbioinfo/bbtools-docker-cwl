package tax;

import java.io.PrintStream;
import java.util.ArrayList;

import server.ServerTools;
import shared.PreParser;
import shared.Timer;
import shared.Tools;

public class TaxClient {
	
	public static void main(String[] args){
		Timer t=new Timer();
		{//Preparse block for help, config files, and outstream
			PreParser pp=new PreParser(args, null, false);
			args=pp.args;
			outstream=pp.outstream;
		}

		ArrayList<Integer> gi=new ArrayList<Integer>();
		ArrayList<String> acc=new ArrayList<String>();
		ArrayList<String> name=new ArrayList<String>();
		ArrayList<String> header=new ArrayList<String>();
		
		boolean slow=true;
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("name") || a.equals("names")){
				if(b!=null){
					for(String s : b.split(",")){
						name.add(s);
					}
				}
			}else if(a.equals("gi")){
				if(b!=null){
					for(String s : b.split(",")){
						gi.add(Integer.parseInt(s));
					}
				}
			}else if(a.equals("accession")){
				if(b!=null){
					for(String s : b.split(",")){
						acc.add(s);
					}
				}
			}else if(a.equals("header")){
				if(b!=null){
					for(String s : b.split(",")){
						header.add(s);
					}
				}
			}else if(a.equals("slow")){
				slow=Tools.parseBoolean(b);
			}else if(a.equals("fast")){
				slow=!Tools.parseBoolean(b);
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		
		if(slow){
			for(String s : name){
				outstream.println(s+"\t"+nameToTaxid(s));
			}
			for(Integer s : gi){
				outstream.println(s+"\t"+giToTaxid(s));
			}
			for(String s : acc){
				outstream.println(s+"\t"+accessionToTaxid(s));
			}
		}else{
			if(name.size()>0){
				int[] ids=nameToTaxidArray(name);
				for(int i=0; i<ids.length; i++){
					outstream.println(name.get(i)+"\t"+ids[i]);
				}
			}
			if(gi.size()>0){
				int[] ids=giToTaxidArray(gi);
				for(int i=0; i<ids.length; i++){
					outstream.println(gi.get(i)+"\t"+ids[i]);
				}
			}
			if(acc.size()>0){
				int[] ids=accessionToTaxidArray(acc);
				for(int i=0; i<ids.length; i++){
					outstream.println(acc.get(i)+"\t"+ids[i]);
				}
			}
		}
		
		t.stopAndPrint();
	}
	
	public static int accessionToTaxid(String accession){
		String s=sendAndRecieve("pt/accession/"+accession);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	public static int giToTaxid(int gi){
		String s=sendAndRecieve("pt/gi/"+gi);
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	public static int nameToTaxid(String name){
		String s=sendAndRecieve("pt/name/"+name.replace(' ', '_'));
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	public static int headerToTaxid(String header){
		String s=sendAndRecieve("pt/name/"+header.replaceAll(" ", "%20"));
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return -1;}
		return Integer.parseInt(s);
	}
	
	
	public static int[] accessionToTaxidArray(String accession){
		String s=sendAndRecieve("pt/accession/"+accession);
		return splitOutput(s);
	}
	
	public static int[] giToTaxidArray(String gi){
		String s=sendAndRecieve("pt/gi/"+gi);
		return splitOutput(s);
	}
	
	public static int[] nameToTaxidArray(String name){
		String s=sendAndRecieve("pt/name/"+name.replace(' ', '_'));
		return splitOutput(s);
	}
	
	public static int[] headerToTaxidArray(String header){
		String s=sendAndRecieve("pt/name/"+header);
		return splitOutput(s);
	}
	
	
	public static int[] accessionToTaxidArray(ArrayList<String> accession){
		String s=sendAndRecieve("pt/accession/"+fuse(accession));
		return splitOutput(s);
	}
	
	public static int[] giToTaxidArray(ArrayList<Integer> gi){
		String s=sendAndRecieve("pt/gi/"+fuse(gi));
		return splitOutput(s);
	}
	
	public static int[] nameToTaxidArray(ArrayList<String> name){
		String s=sendAndRecieve("pt/name/"+fuse(name).replace(' ', '_'));
		return splitOutput(s);
	}
	
	public static int[] headerToTaxidArray(ArrayList<String> header){
		String s=sendAndRecieve("pt/header/"+fuse(header));
		return splitOutput(s);
	}
	
	private static final int[] splitOutput(String s){
		if(s==null || s.length()<1 || !Tools.isDigitOrSign(s.charAt(0))){return null;}
		String[] split=s.split(",");
		int[] ret=new int[split.length];
		for(int i=0; i<split.length; i++){
			ret[i]=Integer.parseInt(split[i]);
		}
		return ret;
	}
	
	private static final String sendAndRecieve(String message){
		final String response;
		if(message.length()<2000){//NERSC Apache limit is around 8kb
			if(message.indexOf(' ')>=0){message=message.replaceAll(" ", "%20");}
			response=ServerTools.sendAndReceive(null, path+message);
		}else{
			response=ServerTools.sendAndReceive(message.getBytes(), path+"$POST");
		}
		return response;
	}
	
	private static String fuse(ArrayList<?> list){
		if(list==null || list.size()<0){return null;}
		StringBuilder sb=new StringBuilder();
		for(Object s : list){
			sb.append(s).append(',');
		}
		sb.setLength(sb.length()-1);
		return sb.toString();
	}
	
	public static String path="https://taxonomy.jgi-psf.org/";
	public static PrintStream outstream=System.err;
	public static boolean verbose=false;
	
}
