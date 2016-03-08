package pacbio;

import java.util.ArrayList;
import java.util.Arrays;

import stream.ConcurrentGenericReadInputStream;
import stream.ConcurrentReadStreamInterface;
import stream.FASTQ;
import stream.FastaReadInputStream;
import stream.Read;


import align2.ListNum;
import align2.Tools;
import dna.Data;
import dna.Timer;
import fileIO.ReadWrite;
import fileIO.FileFormat;
import fileIO.TextStreamWriter;

/**
 * @author Brian Bushnell
 * @date Nov 15, 2012
 *
 */
public class PartitionReads {

	public static void main(String[] args){
		System.err.println("Executing "+(new Object() { }.getClass().getEnclosingClass().getName())+" "+Arrays.toString(args)+"\n");
		
		FastaReadInputStream.SPLIT_READS=false;
		
		Timer t=new Timer();
		t.start();
		
		boolean verbose=false;
		int ziplevel=-1;

		String in1=null;
		String in2=null;
		long maxReads=-1;
		
		String outname1=null;
		String outname2=null;
		FASTQ.PARSE_CUSTOM=false;
		
		for(int i=0; i<args.length; i++){
			final String arg=args[i];
			final String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : "true";
			if("null".equalsIgnoreCase(b)){b=null;}
//			System.err.println("Processing "+args[i]);
			
			if(arg.startsWith("-Xmx") || arg.startsWith("-Xms") || arg.equals("-ea") || arg.equals("-da")){
				//jvm argument; do nothing
			}else if(a.equals("path") || a.equals("root") || a.equals("tempdir")){
				Data.setPath(b);
			}else if(a.equals("fasta") || a.equals("in") || a.equals("input") || a.equals("in1") || a.equals("input1")){
				in1=b;
				if(b.indexOf('#')>-1){
					in1=b.replace("#", "1");
					in2=b.replace("#", "2");
				}
			}else if(a.equals("in2") || a.equals("input2")){
				in2=b;
			}else if(a.startsWith("fastareadlen")){
				FastaReadInputStream.TARGET_READ_LEN=Integer.parseInt(b);
				FastaReadInputStream.SPLIT_READS=(FastaReadInputStream.TARGET_READ_LEN>0);
			}else if(a.startsWith("fastaminread") || a.startsWith("fastaminlen")){
				FastaReadInputStream.MIN_READ_LEN=Integer.parseInt(b);
			}else if(a.equals("fastawrap")){
				FastaReadInputStream.DEFAULT_WRAP=Integer.parseInt(b);
			}else if(a.endsWith("parsecustom")){
				FASTQ.PARSE_CUSTOM=Tools.parseBoolean(b);
				System.out.println("Set FASTQ.PARSE_CUSTOM to "+FASTQ.PARSE_CUSTOM);
			}else if(a.startsWith("partition")){
				partitions=Integer.parseInt(b);
			}else if(a.equals("ziplevel") || a.equals("zl")){
				ReadWrite.ZIPLEVEL=Integer.parseInt(b);
			}else if(a.equals("overwrite") || a.equals("ow")){
				OVERWRITE=Tools.parseBoolean(b);
				System.out.println("Set OVERWRITE to "+OVERWRITE);
			}else if(a.equals("reads") || a.equals("maxreads")){
				maxReads=Long.parseLong(b);
			}else if(a.equals("out") || a.equals("out1")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none") || split.length==1){
					System.out.println("No output file.");
					outname1=null;
				}else{
					outname1=b;
					assert(!outname1.equalsIgnoreCase(outname2));
				}
			}else if(a.equals("out2")){
				if(b==null || b.equalsIgnoreCase("null") || b.equalsIgnoreCase("none") || split.length==1){
					outname2=null;
				}else{
					outname2=b;
					assert(!outname2.equalsIgnoreCase(outname1));
				}
			}else if(a.equals("asciiin") || a.equals("qualityin") || a.equals("qualin") || a.equals("qin")){
				byte ascii_offset=Byte.parseByte(b);
				FASTQ.ASCII_OFFSET=ascii_offset;
				System.out.println("Set fastq input ASCII offset to "+FASTQ.ASCII_OFFSET);
				FASTQ.DETECT_QUALITY=false;
			}else if(a.startsWith("verbose")){
				verbose=Tools.parseBoolean(b);
			}else{
				throw new RuntimeException("Unknown parameter: "+args[i]);
			}
		}
		
		assert(FastaReadInputStream.settingsOK());
		assert(outname1==null || outname1.indexOf('#')>=0 || partitions<2);
		assert(outname2==null || outname2.indexOf('#')>=0 || partitions<2);
		assert(outname1==null || !outname1.equalsIgnoreCase(outname2));
		
		if(in1==null){throw new RuntimeException("Please specify input file.");}
		
		
		final ConcurrentReadStreamInterface cris;
		{
			FileFormat ff1=FileFormat.testInput(in1, FileFormat.FASTQ, null, true, true);
			FileFormat ff2=FileFormat.testInput(in2, FileFormat.FASTQ, null, true, true);
			cris=ConcurrentGenericReadInputStream.getReadInputStream(maxReads, false, true, ff1, ff2);
			if(verbose){System.err.println("Started cris");}
//			Thread th=new Thread(cris);
//			th.start();
		}
		

		TextStreamWriter[] tsw1=new TextStreamWriter[partitions];
		TextStreamWriter[] tsw2=new TextStreamWriter[partitions];
		
		FileFormat ff=FileFormat.testOutput(outname1, FileFormat.FASTQ, null, true, OVERWRITE, false);
		fastq=ff.fastq();
		fasta=ff.fasta();
		bread=ff.bread();
		
		for(int i=0; i<partitions; i++){
			tsw1[i]=new TextStreamWriter(outname1.replaceFirst("#", ""+i), OVERWRITE, false, true);
			if(outname2!=null){
				tsw2[i]=new TextStreamWriter(outname2.replaceFirst("#", ""+i), OVERWRITE, false, true);
			}
		}
		
		long reads=process(tsw1, tsw2, cris);
		t.stop();
		System.out.println("Reads: \t"+reads);
		System.out.println("Time:  \t"+t);
	}
	
	public static long process(TextStreamWriter[] tsw1, TextStreamWriter[] tsw2, ConcurrentReadStreamInterface cris){
		for(TextStreamWriter tsw : tsw1){if(tsw!=null){tsw.start();}}
		for(TextStreamWriter tsw : tsw2){if(tsw!=null){tsw.start();}}
		Thread cristhread=new Thread(cris);
		cristhread.start();
		
		ListNum<Read> ln=cris.nextList();
		ArrayList<Read> readlist=ln.list;
		
		final boolean paired=cris.paired();
		
		long x=0;
		final int div=tsw1.length;
		while(!readlist.isEmpty()){
			
			//System.err.println("Got a list of size "+readlist.size());
			for(int i=0; i<readlist.size(); i++){
				Read r=readlist.get(i);
				if(r!=null){
					final Read r2=r.mate;
					final int mod=(int)(x%div);
					
					StringBuilder a=null, b=null;
					
					if(fastq){
						a=r.toFastq();
						if(paired){b=r2.toFastq();}
					}else if(fasta){
						a=r.toFasta();
						if(paired){b=r2.toFasta();}
					}else if(bread){
						a=r.toText(true);
						if(paired){b=(r2==null ? new StringBuilder(".") : r2.toText(true));}
					}else{
						throw new RuntimeException("Unsupported output format.");
					}
					
					a.append('\n');
					tsw1[mod].print(a);
					if(paired){
						b.append('\n');
						if(tsw2[i]!=null){tsw2[i].print(b);}
						else{tsw1[i].print(b);}
					}
					
					x++;
				}
			}
			
			cris.returnList(ln, readlist.isEmpty());
			
			//System.err.println("Waiting on a list...");
			ln=cris.nextList();
			readlist=ln.list;
		}
		
		//System.err.println("Returning a list... (final)");
		assert(readlist.isEmpty());
		cris.returnList(ln, readlist.isEmpty());
		
		
		for(TextStreamWriter tsw : tsw1){if(tsw!=null){tsw.poison();}}
		for(TextStreamWriter tsw : tsw2){if(tsw!=null){tsw.poison();}}
		ReadWrite.closeStream(cris);
		return x;
	}
	
	public static boolean OVERWRITE=false;
	public static int partitions=2;
	public static boolean fastq=false;
	public static boolean fasta=false;
	public static boolean bread=false;
	
}