import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

public class NewPlayer implements PokerSquaresPlayer{
	
	private final int SIZE = 5; // number of rows/columns in square grid
	private final int NUM_OF_POSITION = SIZE * SIZE; // number of positions in square grid
	private final int NUM_CARDS = Card.NUM_CARDS; // number of cards in deck
	private Random random = new Random(); // pseudorandom number generator for simulation
	private int[] plays = new int[NUM_OF_POSITION]; 
	private int numberofPlays = 0; // number of Cards played into the grid so far
	private PokerSquaresPointSystem system; // point system
	public Card[][] grid = new Card[SIZE][SIZE]; // grid with Card objects or null (for empty positions)
	private Card[] simulationDeck = Card.getAllCards(); // a list of all Cards. As we learn the index of cards in the play deck,
    													// we swap each dealt card to its correct index.  Thus, from index numPlays 
	 													// onward, we maintain a list of undealt cards for simulation.
	private int[][] legalPlayLists = new int[NUM_OF_POSITION][NUM_OF_POSITION]; // stores legal play lists indexed by numPlays (depth)
	
	private final int SIZE_OF_HAND = 5;
	public int StraightFlushScore;
	public int StraightScore;
	public int FlushScore;
	public int TwoPairScore;
	public int ThreePairScore;
	public int FullHouseScore;
	public int RoyalFlushScore;
	public int FourPairScore;
	public int PairScore;
	public int score;
	int cscore=0;
	
	
	@Override
	public void setPointSystem(PokerSquaresPointSystem system, long millis) {
		this.system = system;
		RoyalFlushScore = system.getHandScore(PokerHand.ROYAL_FLUSH);
		StraightFlushScore = system.getHandScore(PokerHand.STRAIGHT_FLUSH);
		FourPairScore = system.getHandScore(PokerHand.FOUR_OF_A_KIND);
		StraightScore = system.getHandScore(PokerHand.STRAIGHT);
		FlushScore = system.getHandScore(PokerHand.FLUSH);
		FullHouseScore = system.getHandScore(PokerHand.FULL_HOUSE);
		ThreePairScore = system.getHandScore(PokerHand.THREE_OF_A_KIND);
		TwoPairScore = system.getHandScore(PokerHand.TWO_PAIR);
		PairScore = system.getHandScore(PokerHand.ONE_PAIR);
	}

	@Override
	public void init() {
		// clear grid
		for (int row = 0; row < SIZE; row++)
			for (int column = 0; column < SIZE; column++)
				grid[row][column] = null;
		// reset numberofPlays
		numberofPlays = 0;
		// (re)initialize list of play positions (row-major ordering)
		for (int i = 0; i < NUM_OF_POSITION; i++)
			plays[i] = i;
		
	}

	@Override
	public int[] getPlay(Card card, long millisRemaining) {
		
		/*
		 * With this algorithm, the player chooses the legal play that has the highest expected score outcome.
		 * This outcome is estimated as follows:
		 *   For each move, many simulated greedy plays to the set depthLimit are performed and the (sometimes
		 *     partially-filled) grid is scored.
		 *   For each greedy play simulation, random undrawn cards are drawn in simulation and the greedy player
		 *     picks a play position that maximizes the score (breaking ties randomly).
		 *   After many such plays, the average score per simulated play is computed.  The play with the highest 
		 *     average score is chosen (breaking ties randomly).   
		 */
		
		// match simDeck to actual play event; in this way, all indices forward from the card contain a list of 
		//   undealt Cards in some permutation.
		
		int indexofcard = numberofPlays;
		while (!card.equals(simulationDeck[indexofcard]))
			indexofcard++;
		simulationDeck[indexofcard] = simulationDeck[numberofPlays];
		simulationDeck[numberofPlays] = card;
		if (numberofPlays < 24) { // not the forced last play
			// compute average time per move evaluation
			int remainingPlays = NUM_OF_POSITION - numberofPlays; // ignores triviality of last play to keep a conservative margin for game completion
			long millisPerPlay = millisRemaining / remainingPlays; // dividing time evenly with future getPlay() calls
			long millisPerMoveEval = millisPerPlay / remainingPlays; // dividing time evenly across moves now considered
			// copy the play positions (row-major indices) that are empty
			System.arraycopy(plays, numberofPlays, legalPlayLists[numberofPlays], 0, remainingPlays);
			double maxAverageScore = Double.NEGATIVE_INFINITY; // maximum average score found for moves so far
			ArrayList<Integer> bestPlays = new ArrayList<Integer>();  // all plays yielding the maximum average score
			for (int i = 0; i < remainingPlays; i++) { // for each legal play position
				int play = legalPlayLists[numberofPlays][i];
				long startTime = System.currentTimeMillis();
				long endTime = startTime + millisPerMoveEval; // compute when MC simulations should end
				makePlay(card, play / SIZE, play % SIZE);  // play the card at the empty position
				int simulationCount = 0;
				int scoreofTotal = 0;
				while (System.currentTimeMillis() < endTime) { // perform as many MC simulations as possible through the allotted time
					// Perform a Monte Carlo simulation of greedy play to the depth limit or game end, whichever comes first.
					try {
					scoreofTotal += simGreedyPlay(32); // accumulate MC simulation scores
					}
					catch(Exception e){
						e.printStackTrace();
					}
					simulationCount++; // increment count of MC simulations
				}
				undoPlay(); // undo the play under evaluation
				// update (if necessary) the maximum average score and the list of best plays
				double averageScore = (double) scoreofTotal / simulationCount;
				if (averageScore >= maxAverageScore) {
					if (averageScore > maxAverageScore)
						bestPlays.clear();
					bestPlays.add(play);
					maxAverageScore = averageScore;
				}
			}
			int bestPlay = bestPlays.get(random.nextInt(bestPlays.size())); // choose a best play (breaking ties randomly)
			// update our list of plays, recording the chosen play in its sequential position; all onward from numPlays are empty positions
			int playIndex = numberofPlays;
			while (plays[playIndex] != bestPlay)
				playIndex++;
			plays[playIndex] = plays[numberofPlays];
			plays[numberofPlays] = bestPlay;
		}
		int[] playPosition = {plays[numberofPlays] / SIZE, plays[numberofPlays] % SIZE}; // decode it into row and column
		makePlay(card, playPosition[0], playPosition[1]); // make the chosen play (not undoing this time)
		return playPosition; // return the chosen play
		
	}
	
	
	/**
	 * From the chosen play, perform simulated Card draws and greedy placement (depthLimit) iterations forward 
	 * and return the resulting grid score.
	 * @param depthLimit - how many simulated greedy plays to perform
	 * @return resulting grid score after greedy MC simulation to given depthLimit
	 */
	private int simGreedyPlay(int depthLimit) {
		if (depthLimit == 0) { // with zero depth limit, return current score
			return system.getScore(grid);
		}
		else { // up to the non-zero depth limit or to game end, iteratively make the given number of greedy plays 0
			int score = Integer.MIN_VALUE;
			int maxScore = Integer.MIN_VALUE;
			int depth = Math.min(depthLimit, NUM_OF_POSITION - numberofPlays); // compute real depth limit, taking into account game end
			for (int d = 0; d < depth; d++) {
				// generate a random card draw
				int ca = random.nextInt(NUM_CARDS - numberofPlays) + numberofPlays;
				Card card = simulationDeck[ca];
				// iterate through legal plays and choose the best greedy play
				int remainingPlays = NUM_OF_POSITION - numberofPlays;
				System.arraycopy(plays, numberofPlays, legalPlayLists[numberofPlays], 0, remainingPlays);
				maxScore = Integer.MIN_VALUE;
				ArrayList<Integer> bestPlays = new ArrayList<Integer>();
				for (int i = 0; i < remainingPlays; i++) {
					int play = legalPlayLists[numberofPlays][i];
					if(numberofPlays>1)
					{
						makePlay(card, play / SIZE, play % SIZE);
						score= Checker(grid);
					}
					else {
					makePlay(card, play / SIZE, play % SIZE);
					}
					
					if (score >= maxScore) {
						if (score > maxScore)
							bestPlays.clear();
						bestPlays.add(play);
						maxScore = score;
					}
					undoPlay();
				}
				int bestPlay = bestPlays.get(random.nextInt(bestPlays.size()));
				makePlay(card, bestPlay / SIZE, bestPlay % SIZE);
			}
			// At this point, the last maxScore value is the end value of this Monte Carlo situation.
						// Undo MC plays.
			for (int d = 0; d < depth; d++) {
				undoPlay();
			}
			return maxScore;
		}
	}

	@Override
	public String getName() {
		return "NewPlayer";
	}
	
	public void makePlay(Card card, int row, int col) {
		// match simDeck to event
		int cardIndex = numberofPlays;
		while (!card.equals(simulationDeck[cardIndex]))
			cardIndex++;
		simulationDeck[cardIndex] = simulationDeck[numberofPlays];
		simulationDeck[numberofPlays] = card;
		
		// update plays to reflect chosen play in sequence
		grid[row][col] = card;
		int play = row * SIZE + col;
		int k = 0;
		while (plays[k] != play)
			k++;
		plays[k] = plays[numberofPlays];
		plays[numberofPlays] = play;
		
		// increment the number of plays taken
		numberofPlays++;
	}

	public void undoPlay() { // undo the previous play
		numberofPlays--;
		int play = plays[numberofPlays];
		grid[play / SIZE][play % SIZE] = null;	
	}

	public int Checker(Card [][]grid) {
		Card rowhand[]= new Card [SIZE_OF_HAND];
		Card colhand[]= new Card [SIZE_OF_HAND];
		int totscore = 0;
		
		// calculating score for row hand
		for (int row = 0; row < SIZE; row++) {
			for (int col = 0; col < SIZE; col++) {
				rowhand[col] = grid[row][col];
			}
			
			if(isStraightFlush(rowhand))
				totscore = totscore + StraightFlushScore;
			if(isStraight(rowhand)) 
				totscore = totscore + StraightScore;
			if(isFlush(rowhand)) 
				totscore = totscore + FlushScore;
			if(isFourPair(rowhand))
				totscore = totscore + FourPairScore;
			if(isFullHouse(rowhand))
				totscore = totscore + FullHouseScore;
			if(isThreePair(rowhand))
				totscore = totscore + ThreePairScore;
			if(isTwoPair(rowhand)) 
				totscore = totscore + TwoPairScore;
			
		}
		
		// calculating score for column hand
		for (int col = 0; col < SIZE; col++) {
			for (int row = 0; row < SIZE; row++) {
				colhand[row] = grid[row][col];
			}
			if(isStraightFlush(colhand))
				totscore = totscore + StraightFlushScore;
			if(isStraight(colhand)) 
				totscore = totscore + StraightScore;
			if(isFlush(colhand)) 
				totscore = totscore + FlushScore;
			if(isFourPair(colhand))
				totscore = totscore + FourPairScore;
			if(isFullHouse(colhand))
				totscore = totscore + FullHouseScore;
			if(isThreePair(colhand))
				totscore = totscore + ThreePairScore;
			if(isTwoPair(colhand)) 
				totscore = totscore + TwoPairScore;
		}
		return totscore;
	}
	
	public static boolean isStraightFlush(Card[] hand)
	{
		if (isFlush(hand) && isStraight(hand))
			return true;
		else
			return false;
	}
	public static boolean isFlush(Card[] hand)
	{
		boolean hasFlush=false;
		int suitC[] = new int [Card.NUM_SUITS];
		for (Card card:hand)
			if (card != null)
                suitC[card.getSuit()]++;
        for(int i =0;i<Card.NUM_SUITS;i++)
        	if (suitC[i] != 0) 
				if (suitC[i] == hand.length)
					hasFlush=true;
        return (hasFlush);
	}
	
	public static boolean isStraight(Card[] hand)
	{
		boolean a1;
		int []ranksC = new int[Card.NUM_RANKS];
        for (Card card:hand) 
            if (card != null) 
                ranksC[card.getRank()]++;
        
        int rank = 0;
		while (rank <= Card.NUM_RANKS - 5 && ranksC[rank] == 0)
			rank++;
		a1 = (rank <= Card.NUM_RANKS - 5 && ranksC[rank] == 1 && ranksC[rank + 1] == 1 && ranksC[rank + 2] == 1 && ranksC[rank + 3] == 1 && ranksC[rank + 4] == 1);
		return a1;
	}
	
	public static boolean isFullHouse(Card[] hand)
	{
		boolean a1;
		int[] ranks = new int[5];
		for (int i=0;i<5;i++)
		{
			if(hand[i]==null)
				ranks[i] = 0;
			else
				ranks[i] = hand[i].getRank();
		}
		Arrays.sort(ranks);
		a1 = ranks[0] == ranks[1]&&
			 ranks[1] == ranks[2]&&
			 ranks[2] == ranks[3]&&
			 ranks[3] == ranks[4];
		
		return(a1);
	}
	public static boolean isFourPair(Card[] hand)
	{
		boolean a1,a2;
		int[]ranks = new int[5];
		for (int i=0; i<5;i++)
		{
			if(hand[i] == null)
				ranks[i] = 0;
			else
				ranks[i] = hand[i].getRank();
		}
		Arrays.sort(ranks);
		a1 = ranks[0] == ranks[1]&&
			 ranks[1] == ranks[2]&&
			 ranks[2] == ranks[3];
		a2 = ranks[1] == ranks[2]&&
			 ranks[2] == ranks[3]&&
			 ranks[3] == ranks[4];
		return (a1||a2);
	}
	public static boolean isThreePair(Card[] hand)
	{
		boolean a1,a2,a3;
		int []ranks = new int[5];
		if ( isFourPair(hand) || isFullHouse(hand) )
	         return(false);
		for (int i =0;i<5;i++)
		{
			if(hand[i] == null)
				ranks[i]=0;
			else
				ranks[i]=hand[i].getRank();
		}
		Arrays.sort(ranks);
		a1 = ranks[0] == ranks[1]&&
			 ranks[1] == ranks[2];
		a2 = ranks[1] == ranks[2]&&
			 ranks[2] == ranks[3];
		a3 = ranks[2] == ranks[3]&&
			 ranks[3] == ranks[4];
		return(a1||a2||a3);
	}
	
	public static boolean isTwoPair(Card[] hand)
	   {
	      boolean a1, a2, a3;
	     if ( isFourPair(hand) || isFullHouse(hand) || isThreePair(hand) )
	         return(false);           
	      int[] ranks = new int[5];
	        for (int i = 0; i < 5; i++) {
	            if (hand[i] == null) {
	                ranks[i] = 0;
	            } else {
	                ranks[i] = hand[i].getRank();
	            }
	        }
	        Arrays.sort(ranks);
                     
	      a1 = ranks[0] == ranks[1] &&
	    	   ranks[2] == ranks[3] ;
	      a2 = ranks[0] == ranks[1] &&
	           ranks[3] == ranks[4] ;
	      a3 = ranks[1] == ranks[2] &&
	           ranks[3] == ranks[4] ;

	      return( a1 || a2 || a3 );
	   }
	
	/*
	public static boolean isPair(Card[] hand)
	{
		boolean a1,a2,a3,a4;
		if ( isFourPair(hand) || isFullHouse(hand) || isThreePair(hand) || isTwoPair(hand) )
	         return(false);           
	      int[] ranks = new int[5];
	        for (int i = 0; i < 5; i++) {
	            if (hand[i] == null) {
	                ranks[i] = 0;
	            } else {
	                ranks[i] = hand[i].getRank();
	            }
	        }
	        Arrays.sort(ranks);
	        a1 = ranks[0] == ranks[1];
	 	    a2 = ranks[1] == ranks[2];
	 	    a3 = ranks[2] == ranks[3];
	 	    a4 = ranks[3] == ranks[4];
	 	  return(a1||a2||a3||a4);
	}
	*/
	
	public static void main(String[] args) {
		
		PokerSquaresPointSystem system = PokerSquaresPointSystem.getBritishPointSystem();
		System.out.println(system);
		new PokerSquares(new NewPlayer(), system).play(); // play a single game
	}

}