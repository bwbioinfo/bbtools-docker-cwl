package assemble;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicIntegerArray;

import dna.AminoAcid;
import kmer.AbstractKmerTableSet;
import shared.Tools;
import structures.ByteBuilder;
import ukmer.AbstractKmerTableU;
import ukmer.HashArrayU1D;
import ukmer.HashForestU;
import ukmer.Kmer;
import ukmer.KmerNodeU;
import ukmer.KmerTableSetU;

/**
 * Designed for removal of dead ends (aka hairs).
 * @author Brian Bushnell
 * @date Jun 26, 2015
 *
 */
public class Shaver2 extends Shaver {
	
	
	/*--------------------------------------------------------------*/
	/*----------------         Constructor          ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public Shaver2(KmerTableSetU tables_, int threads_){
		this(tables_, threads_, 1, 1, 1, 1, 3, 100, 100, true, true);
	}
	
	public Shaver2(KmerTableSetU tables_, int threads_,
			int minCount_, int maxCount_, int minSeed_, int minCountExtend_, float branchMult2_, int maxLengthToDiscard_, int maxDistanceToExplore_,
			boolean removeHair_, boolean removeBubbles_){
		super(tables_, threads_, minCount_, maxCount_, minSeed_, minCountExtend_, branchMult2_, maxLengthToDiscard_, maxDistanceToExplore_, removeHair_, removeBubbles_);
		tables=tables_;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/

	@Override
	final AbstractExploreThread makeExploreThread(int id_){return new ExploreThread(id_);}
	@Override
	final AbstractShaveThread makeShaveThread(int id_){return new ShaveThread(id_);}
	
	
	/*--------------------------------------------------------------*/
	/*----------------       Dead-End Removal       ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public boolean exploreAndMark(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount,
			int maxLengthToDiscard, int maxDistanceToExplore, boolean prune
			, long[][] countMatrixT, long[][] removeMatrixT
			){
		bb.clear();
		assert(kmer.len>=kmer.kbig);
		if(findOwner(kmer)>STATUS_UNEXPLORED){return false;}
		
		bb.appendKmer(kmer);
//		assert(tables.getCount(kmer)==1) : tables.getCount(kmer);//count>0 && count<=maxCount
		final int a=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);
		
		bb.reverseComplementInPlace();
		kmer=tables.rightmostKmer(bb, kmer);
//		assert(tables.getCount(kmer)==1) : tables.getCount(kmer);//count>0 && count<=maxCount
		final int b=explore(kmer, bb, leftCounts, rightCounts, minCount, maxCount, maxDistanceToExplore);

		final int min=Tools.min(a, b);
		final int max=Tools.max(a, b);
		
		countMatrixT[min][max]++;
		
		if(a==TOO_LONG || a==TOO_DEEP || a==LOOP || a==FORWARD_BRANCH){
			claim(bb, STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(b==TOO_LONG || b==TOO_DEEP || b==LOOP || b==FORWARD_BRANCH){
			claim(bb, STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(bb.length()-kbig>maxLengthToDiscard){
			claim(bb, STATUS_EXPLORED, false, kmer);
			return false;
		}
		
		if(removeHair && min==DEAD_END){
			if(max==DEAD_END || max==BACKWARD_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false, kmer);
				if(verbose || verbose2){System.err.println("Claiming ("+a+","+b+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		if(removeBubbles){
			if(a==BACKWARD_BRANCH && b==BACKWARD_BRANCH){
				removeMatrixT[min][max]++;
				boolean success=claim(bb, STATUS_REMOVE, false, kmer);
				if(verbose || verbose2){System.err.println("Claiming ("+a+","+b+") length "+bb.length()+": "+bb);}
				assert(success);
				return true;
			}
		}
		
		claim(bb, STATUS_EXPLORED, false, kmer);
		return false;
	}
	
	/** Explores a single unbranching path in the forward direction.
	 * Returns reason for ending in this direction:
	 *  DEAD_END, TOO_LONG, TOO_DEEP, FORWARD_BRANCH, BACKWARD_BRANCH */
	public int explore(Kmer kmer, ByteBuilder bb, int[] leftCounts, int[] rightCounts, int minCount, int maxCount, int maxLength0){
		if(verbose){outstream.println("Entering explore with bb.length()="+bb.length());}
		assert(bb.length()>=kmer.kbig && kmer.len>=kmer.kbig);
		if(bb.length()==0){bb.appendKmer(kmer);}
		
		final int initialLength=bb.length();
		final int maxLength=maxLength0+kbig;
		
		final long firstKey=kmer.xor();
		HashArrayU1D table=tables.getTable(kmer);
		int count=table.getValue(kmer);
//		kmer.verify(false);
//		kmer.verify(true);
		assert(count>=minCount && count<=maxCount) : count+", "+Kmer.MASK_CORE+", "+kmer.verify(false)+", "+kmer.verify(true);
		
		int nextRightMaxPos=fillRightCounts(kmer, rightCounts);
		int nextRightMax=rightCounts[nextRightMaxPos];
		if(nextRightMax<minCount){
			if(verbose){outstream.println("Returning DEAD_END: rightMax="+nextRightMax);}
			return DEAD_END;
		}
		
		while(bb.length()<=maxLength){
			
			final int rightMaxPos=nextRightMaxPos;
			final int rightMax=rightCounts[rightMaxPos];
			final int rightSecondPos=Tools.secondHighestPosition(rightCounts);
			final int rightSecond=rightCounts[rightSecondPos];
			
			if(verbose){
				outstream.println("kmer: "+toText(kmer));
				outstream.println("Right counts: "+count+", "+Arrays.toString(rightCounts));
				outstream.println("rightMaxPos="+rightMaxPos);
				outstream.println("rightMax="+rightMax);
				outstream.println("rightSecondPos="+rightSecondPos);
				outstream.println("rightSecond="+rightSecond);
			}
			
			final int prevCount=count;
			
			//Generate the new base
			final byte b=AminoAcid.numberToBase[rightMaxPos];
			final long x=rightMaxPos;
			long evicted=kmer.addRightNumeric(x);
			
			//Now consider the next kmer
			if(kmer.xor()==firstKey){
				if(verbose){outstream.println("Returning LOOP");}
				return LOOP;
			}
			table=tables.getTable(kmer);
			
			assert(table.getValue(kmer)==rightMax);
			count=rightMax;
			
			
			{//Fill right and look for dead end
				nextRightMaxPos=fillRightCounts(kmer, rightCounts);
				nextRightMax=rightCounts[nextRightMaxPos];
				if(nextRightMax<minCount){
					if(verbose){outstream.println("Returning DEAD_END: rightMax="+rightMax);}
					return DEAD_END;
				}
			}
			
			
			{//Look left
				final int leftMaxPos=fillLeftCounts(kmer, leftCounts);
				final int leftMax=leftCounts[leftMaxPos];
				final int leftSecondPos=Tools.secondHighestPosition(leftCounts);
				final int leftSecond=leftCounts[leftSecondPos];
				
//				assert(leftMax==1 || leftMax==0) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
				
				if(verbose){
					outstream.println("Left counts: "+count+", "+Arrays.toString(leftCounts));
					outstream.println("leftMaxPos="+leftMaxPos);
					outstream.println("leftMax="+leftMax);
					outstream.println("leftSecondPos="+leftSecondPos);
					outstream.println("leftSecond="+leftSecond);
				}
				
				if(leftSecond>=minCount || leftMax>prevCount){//Backward branch
//					assert(leftSecond==1) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
					if(leftMax>prevCount){
						if(verbose){outstream.println("Returning BACKWARD_BRANCH_LOWER: " +
								"count="+count+", prevCount="+prevCount+", leftMax="+leftMax+", leftSecond="+leftSecond);}
						return BACKWARD_BRANCH;
					}else{
						assert(leftMax==prevCount);
						if(leftMax>=2*leftSecond){//This constant is adjustable
//							assert(false) : prevCount+" -> "+Arrays.toString(leftCounts)+", "+count+", "+Arrays.toString(rightCounts);
							//keep going
						}else{
							if(verbose){outstream.println("Returning BACKWARD_BRANCH_SIMILAR: " +
								"count="+count+", prevCount="+prevCount+", leftMax="+leftMax+", leftSecond="+leftSecond);}
							return BACKWARD_BRANCH;
						}
					}
				}
				
			}
			
			//Look right
			if(rightSecond>=minCount){
				if(verbose){outstream.println("Returning FORWARD_BRANCH: rightSecond="+rightSecond);}
				return FORWARD_BRANCH;
			}
			
			if(count>maxCount){
				if(verbose){outstream.println("Returning TOO_DEEP: rightMax="+rightMax);}
				return TOO_DEEP;
			}
			
			bb.append(b);
			if(verbose){outstream.println("Added base "+(char)b);}
		}
		
		assert(bb.length()>maxLength);
		if(verbose){outstream.println("Returning TOO_LONG: length="+bb.length());}
		return TOO_LONG;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         ExploreThread        ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Searches for dead ends.
	 */
	class ExploreThread extends AbstractExploreThread{
		
		/**
		 * Constructor
		 */
		public ExploreThread(int id_){
			super(id_, kbig);
		}
		
		@Override
		boolean processNextTable(final Kmer kmer, Kmer temp){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final int[] counts=table.values();
			final int max=counts.length;
//			final int max=table.arrayLength();
			if(startFromHighCounts){
//				for(int cell=0; cell<max; cell++){
//					int x=processCell_high(table, cell, kmer, temp);
//					deadEndsFoundT+=x;
//				}
				for(int cell=0; cell<max; cell++){//Not noticeably faster
					final int count=counts[cell];
					if(count>maxCount){
						int x=processCell_high(table, cell, kmer, temp, count);
						deadEndsFoundT+=x;
					}
				}
			}else{
				for(int cell=0; cell<max; cell++){
					int x=processCell_low(table, cell, kmer, temp);
					deadEndsFoundT+=x;
				}
			}
			return true;
		}
		
		@Override
		boolean processNextVictims(final Kmer kmer, Kmer temp){
			final int tnum=nextVictims.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
			final HashArrayU1D table=tables.getTable(tnum);
			final HashForestU forest=table.victims();
			final int max=forest.arrayLength();
			if(startFromHighCounts){
				for(int cell=0; cell<max; cell++){
					KmerNodeU kn=forest.getNode(cell);
					int x=traverseKmerNodeU_high(kn, kmer, temp);
					deadEndsFoundT+=x;
				}
			}else{
				for(int cell=0; cell<max; cell++){
					KmerNodeU kn=forest.getNode(cell);
					int x=traverseKmerNodeU_low(kn, kmer, temp);
					deadEndsFoundT+=x;
				}
			}
			return true;
		}
		
		private int traverseKmerNodeU_high(KmerNodeU kn, Kmer kmer, Kmer temp){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNodeU_high(kn, kmer, temp);
				if(kn.left()!=null){
					sum+=traverseKmerNodeU_high(kn.left(), kmer, temp);
				}
				if(kn.right()!=null){
					sum+=traverseKmerNodeU_high(kn.right(), kmer, temp);
				}
			}
			return sum;
		}
		
		private int traverseKmerNodeU_low(KmerNodeU kn, Kmer kmer, Kmer temp){
			int sum=0;
			if(kn!=null){
				sum+=processKmerNodeU_low(kn, kmer, temp);
				if(kn.left()!=null){
					sum+=traverseKmerNodeU_low(kn.left(), kmer, temp);
				}
				if(kn.right()!=null){
					sum+=traverseKmerNodeU_low(kn.right(), kmer, temp);
				}
			}
			return sum;
		}
		
		/*--------------------------------------------------------------*/
		
		//old
		private int processCell_low(HashArrayU1D table, int cell, Kmer kmer0, Kmer temp){
			int count=table.readCellValue(cell);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=table.getCellOwner(cell);
			if(owner>STATUS_UNEXPLORED){return 0;}
			Kmer kmer=table.fillKmer(cell, kmer0);
			if(kmer==null){return 0;}
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+kmer);}
			
			return processKmer_low(kmer, temp);
		}
		
		//old
		private int processKmerNodeU_low(KmerNodeU kn, Kmer kmer, Kmer temp){
			kmer.setFrom(kn.pivot());
			final int count=kn.getValue(kmer);
			if(count<minSeed || count>maxCount){return 0;}
			int owner=kn.getOwner(kmer);
			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer_low(kmer, temp);
		}
		
		//old
		private int processKmer_low(Kmer original, Kmer temp){
			kmersTestedT++;
			boolean b=exploreAndMark(original, builderT, leftCounts, rightCounts, minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true
					, countMatrixT, removeMatrixT
				);
			return b ? 1 : 0;
		}
		

		
		/*--------------------------------------------------------------*/
		
		//new
		private int processCell_high(HashArrayU1D table, int cell, Kmer kmer0, Kmer temp, int count){
			
//		private int processCell_high(HashArrayU1D table, int cell, Kmer kmer0, Kmer temp){
//			int count=table.readCellValue(cell);
			
			if(count<=maxCount){return 0;}
//			if(shaveFast && maxCount==1 && count>=6){return 0;}
//			int owner=table.getCellOwner(cell);
//			if(owner>STATUS_UNEXPLORED){return 0;}
			Kmer kmer=table.fillKmer(cell, kmer0);
			if(kmer==null){return 0;}
			if(verbose){outstream.println("id="+id+" processing cell "+cell+"; \tkmer="+kmer);}
//			assert(false) : shaveFast;
			
			return processKmer_high(kmer, temp, count);
		}
		
		//new
		private int processKmerNodeU_high(KmerNodeU kn, Kmer kmer, Kmer temp){
			kmer.setFrom(kn.pivot());
			final int count=kn.getValue(kmer);
			if(count<=maxCount){return 0;}
//			if(shaveFast && maxCount==1 && count>=6){return 0;}
//			int owner=kn.getOwner(kmer);
//			if(owner>STATUS_UNEXPLORED){return 0;}
			
			return processKmer_high(kmer, temp, count);
		}
		
		//new
		private int processKmer_high(final Kmer original, Kmer kmer, final int count0){
			int sum=0;
			
//			assert(false) : "This loop is broken and removes kmers with counts outside the range, including non-existing kmers.";
//			for(long i=0; i<4; i++){
//				kmer.setFrom(original);
//				long old=kmer.addRightNumeric(i);
//				HashArrayU1D table=(HashArrayU1D) tables.getTable(kmer);
//				int count=table.getCount(kmer);
//				if(count>0 && count<=maxCount && tables.getOwner(kmer)<=STATUS_UNEXPLORED){
//					kmersTestedT++;
//					boolean b=exploreAndMark(kmer, builderT, leftCounts, rightCounts,
//							minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
//							countMatrixT, removeMatrixT);
//					if(b){sum++;}
//				}
//			}
			
			assert(original.len()==kbig);
			
			//This loop appears to work correctly.
			sum+=processKmer_high_leftLoop(original, kmer, count0);
			
			original.rcomp();
			assert(original.len()==kbig);
			sum+=processKmer_high_leftLoop(original, kmer, count0);
			return sum;
		}
		
		private int processKmer_high_leftLoop(final Kmer original, Kmer kmer, final int count0){
			if(shaveVFast){//Optional accelerator; only examines kmers that actually branch
				int inf=Integer.MAX_VALUE;
				int maxHighBranch=-1, minHighBranch=inf, highBranches=0;
				int maxLowBranch=-1, minLowBranch=inf, lowBranches=0;
				for(long i=0; i<4; i++){
					kmer.setFrom(original);
					long old=kmer.addLeftNumeric(i);
					assert(kmer.len()>=kbig) : kmer.len()+", "+kbig;
					HashArrayU1D table=(HashArrayU1D) tables.getTable(kmer);
					int count=table.getValue(kmer);
					if(count>0){
						if(count<=maxCount){
							minLowBranch=Tools.min(count, minLowBranch);
							maxLowBranch=Tools.max(count, maxLowBranch);
							lowBranches++;
						}else{
							minHighBranch=Tools.min(count, minHighBranch);
							maxHighBranch=Tools.max(count, maxHighBranch);
							highBranches++;
						}
					}
				}

				if(maxLowBranch<0){return 0;} //This is fine
//				if(maxHighBranch<0){return 0;} //Speed increase but quality decrease

				if(highBranches+lowBranches<2){return 0;} //Speed increase (25% for shave) but contiguity decrease (~1%).
//				if(maxLowBranch>0 && maxLowBranch*16<minHighBranch){return 0;} //Speed increase but quality decrease
			}
			
			int sum=0;
			for(long i=0; i<4; i++){
				kmer.setFrom(original);
				long old=kmer.addLeftNumeric(i);
				assert(kmer.len()>=kbig) : kmer.len()+", "+kbig;
				HashArrayU1D table=(HashArrayU1D) tables.getTable(kmer);
				int count=table.getValue(kmer);
				if(count>0 && count<=maxCount && table.getOwner(kmer)<=STATUS_UNEXPLORED){
					kmersTestedT++;
					boolean b=exploreAndMark(kmer, builderT, leftCounts, rightCounts,
							minCount, maxCount, maxLengthToDiscard, maxDistanceToExplore, true,
							countMatrixT, removeMatrixT);
					if(b){sum++;}
				}
			}
			return sum;
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          ShaveThread         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * Removes dead-end kmers.
	 */
	private class ShaveThread extends AbstractShaveThread{

		/**
		 * Constructor
		 */
		public ShaveThread(int id_){
			super(id_);
		}
		
		@Override
		boolean processNextTable(){
			final int tnum=nextTable.getAndAdd(1);
			if(tnum>=tables.ways){return false;}
//			long x=0;
			final HashArrayU1D table=tables.getTable(tnum);
			final AtomicIntegerArray owners=table.owners();
			final int[] values=table.values();
			final int max=table.arrayLength();
			for(int cell=0; cell<max; cell++){
				if(owners.get(cell)==STATUS_REMOVE){
//					x++;
					values[cell]=0;
				}
			}
			for(KmerNodeU kn : table.victims().array()){
				if(kn!=null){traverseKmerNodeU(kn);}
			}
			
			table.clearOwnership();
			kmersRemovedT+=table.regenerate(0);
//			outstream.println(x);
			return true;
		}
		
		private void traverseKmerNodeU(KmerNodeU kn){
			if(kn==null){return;}
			if(kn.owner()==STATUS_REMOVE){kn.set(0);}
			traverseKmerNodeU(kn.left());
			traverseKmerNodeU(kn.right());
		}
		
	}
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Recall Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public boolean claim(final ByteBuilder bb, final int id, final boolean exitEarly, Kmer kmer){
		return claim(bb.array, bb.length(), id, exitEarly, kmer);
	}
	
	public boolean claim(final byte[] bases, final int blen, final int id, boolean exitEarly, final Kmer kmer){
		if(blen<kbig){return false;}
		if(verbose){outstream.println("Thread "+id+" claim start.");}
		int len=0;
		kmer.clear();
		boolean success=true;
		/* Loop through the bases, maintaining a forward and reverse kmer via bitshifts */
		for(int i=0; i<blen && success; i++){
			final byte b=bases[i];
			final long x=AminoAcid.baseToNumber[b];
			kmer.addRight(b);

			if(x<0){len=0;}
			else{len++;}
			assert(len==kmer.len);

			if(len>=kbig){
				assert(countWithinLimits(kmer)) : getCount(kmer)+", "+len+", "+i+", "+blen+"\n"+kmer.toString();
				success=claim(kmer, id/*, rid, i*/);
				success=(success || !exitEarly);
			}
		}
		return success;
	}
	
	final boolean countWithinLimits(Kmer kmer){
		int count=getCount(kmer);
		return count>=minCount && count<=maxCount;
	}
	
	int getCount(Kmer kmer){return tables.getCount(kmer);}
	boolean claim(Kmer kmer, int id){return tables.claim(kmer, id);}
	boolean doubleClaim(ByteBuilder bb, int id/*, long rid*/, Kmer kmer){return tables.doubleClaim(bb, id, kmer/*, rid*/);}
//	boolean claim(ByteBuilder bb, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(bb, id/*, rid*/, earlyExit, kmer);}
//	boolean claim(byte[] array, int len, int id, /*long rid, */boolean earlyExit, Kmer kmer){return tables.claim(array, len, id/*, rid*/, earlyExit, kmer);}
	int findOwner(Kmer kmer){return tables.findOwner(kmer);}
	int findOwner(ByteBuilder bb, int id, Kmer kmer){return tables.findOwner(bb, id, kmer);}
	int findOwner(byte[] array, int len, int id, Kmer kmer){return tables.findOwner(array, len, id, kmer);}
	void release(ByteBuilder bb, int id, Kmer kmer){tables.release(bb, id, kmer);}
	void release(byte[] array, int len, int id, Kmer kmer){tables.release(array, len, id, kmer);}
	int fillRightCounts(Kmer kmer, int[] counts){return tables.fillRightCounts(kmer, counts);}
	int fillLeftCounts(Kmer kmer, int[] counts){return tables.fillLeftCounts(kmer, counts);}
	static StringBuilder toText(Kmer kmer){return AbstractKmerTableU.toText(kmer);}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	@Override
	AbstractKmerTableSet tables(){return tables;}
	
	final KmerTableSetU tables;
	
}
