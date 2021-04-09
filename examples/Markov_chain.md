# Gambler's Ruin

Gambler's ruin is a classic probability problem that can be thought as a **Markov Chain**. The problem usually involves a gambler playing games for money. For each winning game the gambler gains some money. Conversely, the gambler loses some money for a loosing game. Being the gambler, he/she is obsessed with gambling, and the gambler quits only after either winning the entire table or losing all of his/her money. A variation of the Gambler's ruin can look at this:

> The Gambler has $1, and the opponent has $2. For each game the winner gets $1 from the other. As a better player, the gambler wins 2/3 of the games. They would play until either one is bankrupt. What is probability of the gambler winning?

## Approach

The first step is to pick the state space and define the transition probabilities. It is fairly straight forward for this simple problem. **The state can be the gambler's current money count.** The transition probability then looks at:

```text
Transition_Probability =  |1    0    0    0  |
                          |1/3  0    2/3  0  |
                          |0    1/3  0    2/3|
                          |0    0    0    1  |
```

The gambler starts at **State 1**. The next state is **State 0** with probability 1/3, and **State 2** with probability 2/3. Similarly, the gambler can get from **State 2** to **State 1** with probability 1/3, and from **State 2** to **State 3** with probability 2/3. Once the gambler is in **State 3**, the game ends because he/she wins all of the opponent's money. The game can also end when the gambler is in **State 0** as he/she becomes bankrupt.

## Results
Unsurprisingly the gambler, as the better player, is more likely to win the entire table than to lose it all. The probability of winning it all edges out slightly at 4/7 (~0.571). Conversely the probability of bankrupting is 3/7.


### Code

```cpp
/* Gambler's ruin */
// The '1' in diagnal is optional.
// std::string graphStr {"BA1,BC2,CB1,CD2"};

std::string graphStr {"AA1,BA1,BC2,CB1,CD2,DD1"};
Adjmat<int> myAdjmat = convert2DirectedAdjmat<int>(graphStr);

// Normalize the rows to represent probability.
auto normAdjmat = normalizeRow(myAdjmat);
DWGraph myg(normAdjmat);
Markovchain mc(myg);

// Find transition probability from State B to all terminal states.
auto prob_map = mc.find_transit_prob(1);
printf("Probability of winning the table : %g\n", prob_map[3]);
printf("Probability of losing all money : %g\n", prob_map[0]);

```
