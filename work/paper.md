# Lucrehulk
### CCU Inflation Research Paper

---

## 1. Objective & Table of Contents

### 1.1 Objective

The objective of this paper is to design a mathematical system for inflating the displayed concurrent user count (CCU) of a multiplayer lobby with fake players in a way that reduces the probability of detection. This will be based on designing a function $F$ that returns a whole number, which is the number of fake players to add to a lobby. The displayed count $P + F$ should be as large as possible without the extra fake players being noticed by real players. 

For example, assume a simple 2D game with a lobby where a fixed-size, roamable arena spans $X_r \times Y_r$ units as float values. Assume players are circles with an arbitrary radius $R$, each of which has a field of view of dimensions $X_{fov} \times Y_{fov}$. Assume there is a lobby CCU display in each player's view. The question becomes: what kind of model describes how many fake players can be added to the CCU display without real players noticing, and how quickly can these fake players be added without real players noticing? This is the form of question this paper intends to answer.

To answer this question, two different mathematical systems, one for steady state systems that only account for inflation based on players, and a complete system that can handle transient phases and accounts for players and time elapsed, will be developed. To develop these, two separate forms of "detectability" (how easily detectable inflation is) will be analyzed, those being population and temporal.

Additionally, using Rust, a simulation will be developed to help experimentally determine empirical values within the inflation system derived here. Exact experiments, however, will not be designed in this paper. 

This work can be used by researchers to learn more about how bad actors may fraudulently inflate populations.

### 1.2 Table of Contents

| Section | Purpose |
|---|---|
| 2 | Notation |
| 3 | Defining population detectability $D_P$ |
| 4 | Deriving lobby and information awareness |
| 5 | Deriving a model of the form of $F(P)$ |
| 6 | Extending the model to $F(P, t)$ |
| 7 | Final system compilation |
| 8 | Extension to multiple lobbies with an aggregate CCU |
| 9 | Simulation justification for experiments |
| 10 | Temporal simulation design |
| 11 | Summary of findings |

## 2. Notation

### 2.1 Model symbols

| Symbol | Meaning | Range |
|---|---|---|
| $P$ | True player count | integer, $\geq 0$ |
| $F$ | Number of fake players added | integer, $\geq 0$ |
| $\phi$ | Inflated fraction, $F/(P+F)$ | $[0, 1)$ for $P \gt 0$ |
| $t$ | Time since inflation began (in Section 10, the tick index within a trial, so time since inflation began is $t - \text{start}$) | $t \geq 0$ |
| $D_P$ | Population detectability | $[0, 1]$ |
| $D_T$ | Temporal detectability, a rate per tick | $[0, 1]$ |
| $L_A$ | Per-user lobby awareness | $[0, 1]$ |
| $L_K$ | Lobby knowledge | $[0, 1]$ |
| $L_C$ | Lobby chaos | $[0, 1]$ |
| $K_C$ | Accumulated knowledge | $[0, 1]$ |
| $U_A$ | User awareness | $[0, 1]$ |
| $V_A$, $R_A$, $V_F$ | Viewable area, total lobby area, view fraction $V_A / R_A$ | $V_F \in [0, 1]$ |
| $w_{V_F}$, $w_{K_C}$ | Weights on $V_F$ and $K_C$ in $L_K$ | $\geq 0$, not both zero |
| $I_C$ | Information completeness | $[0, 1]$ |
| $I_P$ | Information processability | $[0, 1]$ |
| $h$ | Processability function, $\mathbb{E}[I_P] = \mathbb{E}[h(L_C, U_A)]$, decreasing in $L_C$ and increasing in $U_A$ | $[0, 1]$ |
| $I_A$ | Total information awareness | $[0, 1]$ |
| $w_{I_C}$, $w_{I_P}$ | Weights on $I_C$ and $I_P$ | $\geq 0$, not both zero |
| $k$ | Inflation ratio, $F = kP$ at steady state | $k \gt 0$, and $k \leq 1$ in practice (Section 6.5) |
| $k_{\max}$ | Largest allowed $k$, $\frac{D_{P,0}}{1 - D_{P,0}}$ | $(0, 1]$ in practice |
| $T_0$ | Damping duration | whole number of ticks, $T_0 \geq 1$ |
| $T_{0,\min}$ | Smallest allowed $T_0$, $\frac{k}{D_{T,0}}$ (the smallest whole $T_0$ meeting it is $\lceil T_{0,\min} \rceil$) | $\gt 0$ |
| $D_{P,0}$, $D_{T,0}$ | Detectability limits, found experimentally | $(0, 1)$, with $D_{P,0} \leq \frac{1}{2}$ when $k \leq 1$ |
| $N$ | Number of lobbies | integer, $N \geq 2$ |
| $P_n$, $F_n$ | Real and fake players in lobby $n$ | integer, $\geq 0$ |
| $P_{total}$, $F_{total}$ | Real and fake players summed across all $N$ lobbies | integer, $\geq 0$ |

A tick is one discrete unit of time (for example, one minute). The quantities $t$ and $T_0$ are plain numbers counted in ticks, so $D_T$ is a rate per tick.

### 2.2 Simulation symbols (Section 10)

| Symbol | Code | Meaning | Range |
|---|---|---|---|
| $M$ | `len` | Trial length in ticks, $M = 2 \cdot T_{max}$ | integer |
| $P_0$ | `P_0` | Starting true player count in a trial (the target average CCU) | $P_{\min} \leq P_0 \leq P_{\max}$ |
| $P_{\min}$, $P_{\max}$ | `P_MIN`, `P_MAX` | Lower and upper bounds on the true count | integer |
| $P_f[t]$ | `fluctuations[t]` | Natural fluctuation applied at tick $t$, so $P[t+1] = P[t] + P_f[t+1]$ | integer |
| $\mu$, $\sigma$ | `MU`, `SIG` | Mean and standard deviation of $P_f[t]$ | $\mu \approx 0$, $\sigma \gt 0$ |
| $P_d(x)$ | `norm(x)` | Density of the fluctuation distribution | $\geq 0$ |
| $k$ | `K` | Inflation ratio, fixed to $k_{\max}$ in actual use (not to be confused with $K_C$) | $0 \lt k \leq 1$ |
| $T_{min}$, $T_{max}$ | `MIN_T`, `MAX_T` | Smallest and largest $T_0$ tested (not the model's $T_{0,\min}$) | whole ticks, $1 \leq T_{min} \leq T_{max}$ |
| $T_{jump}$ | `T_JUMP` | Step between consecutive tested $T_0$ values | whole ticks, $\geq 1$ |
| $N_{inflated}$ | Not Applicable | Number of inflated trials, $\left\lfloor \frac{T_{max} - T_{min}}{T_{jump}} \right\rfloor + 1$ | integer |
| $\rho_C$ | `CONTROL_FRACT` | Control trials as a fraction of inflated trials | $\geq 0$ |
| $N_{control}$ | `control_count` | Number of control trials, $\text{round}(\rho_C \cdot N_{inflated})$ | integer |
| $\text{start}$, $\text{end}$ | `inflation_period_start`, `inflation_period_end` | Tick at which the ramp begins, and tick at which it reaches full value ($\text{end} = \text{start} + T_0$) | $0 \leq \text{start} \leq M - T_0 - 1$ |

---

## 3. Population Detectability

### 3.1 Detectability and the inflated fraction

Detectability is how easy it is to tell that a displayed player count is inflated. Population detectability measures this based only on the current number of players $P$. The main factor is the fraction of the total displayed population $P + F$ that is made up of the inflation count $F$,

$$\phi = \frac{F}{P+F}, \qquad P \gt 0$$

As $\phi$ goes up, a user joining the lobby is more likely to notice the inflation. As it goes down, the user is less likely to notice. Or more specifically, $\phi \uparrow \Rightarrow D_P \uparrow$ and $\phi \downarrow \Rightarrow D_P \downarrow$. Other factors matter too (some of which are explained in later sections, along with how their effects are absorbed and implicitly accounted for), but this fraction is generally a reasonable estimate of population detectability. This gives the approximation:

$$D_P \approx \phi = \frac{F}{P+F}$$

### 3.2 $F$ at a fixed inflated fraction

Holding $\phi$ fixed and solving for $F$:

$$\phi = \frac{F}{P+F}$$

$$\frac{1}{\phi} = \frac{P+F}{F} = \frac{P}{F} + 1$$

$$\frac{1}{\phi} - 1 = \frac{P}{F}$$

$$F = \frac{P}{\frac{1}{\phi}-1} = \frac{\phi P}{1-\phi}$$

So at a given $\phi$, $F$ is directly proportional to $P$. That is, $F \propto P$. A bigger lobby can hold more fake players at the same level of detectability, so

$$P \uparrow \Rightarrow F \uparrow$$

The concept of lobby awareness can now be analyzed to explain more concretely the factors that drive detectability.

---

## 4. Lobby Awareness

### 4.1 Definition

Let $L_A \in [0, 1]$ be how sensitive an individual user is to inflation in a lobby. $L_A = 1$ means the user notices even a single fake player. $L_A = 0$ means the user would never notice inflation, no matter how many players are fake. It is effectively positively related to both forms of detectability, meaning $L_A \uparrow \Rightarrow D_P \uparrow$ and $L_A \uparrow \Rightarrow D_T \uparrow$. Note, though, that the weight of these factors can differ between $D_P$ and $D_T$. For example, factors like lobby knowledge are much less significant in temporal detection analysis. The components below are therefore described mainly in relation to $D_P$, though they still apply to $D_T$ to a lesser degree. Since $L_A$ differs from user to user, its expected value $\mathbb{E}[L_A]$ is used.

### 4.2 Components of lobby awareness

Lobby awareness depends on three major factors:

**Lobby knowledge ($L_K$).** This is how much knowledge of the lobby a user holds at a given point. More specifically, it is the ratio of the user's own knowledge to the entire knowledge available within a lobby. While knowledge is difficult to quantify, comparing a user's field of view (the entire viewable area) to the entire room area may yield a general estimate of what fraction of the total knowledge the user holds. With viewable area $V_A$ and total lobby (room) area $R_A$, the view fraction is $V_F = \frac{V_A}{R_A} \in [0, 1]$.

Additionally, a user retains knowledge from points in the past. This means an individual user's ability to detect inflation increases the longer the user remains in a lobby and accumulates knowledge, or more specifically, $K_C \uparrow \Rightarrow L_K \uparrow \Rightarrow L_A \uparrow$. To account for this, accumulated knowledge $K_C \in [0, 1]$ is introduced as the fraction of the lobby's total information that the user has gathered from earlier views. Total $L_K$ is approximated by a weighted average of $V_F$ and $K_C$ with weights $w_{V_F}, w_{K_C} \geq 0$, not both zero:

$$L_K \approx \frac{w_{V_F} V_F + w_{K_C} K_C}{w_{V_F} + w_{K_C}}$$

Of course, different parts of a lobby hold different amounts of information. On average, though, it is reasonable to assume information is spread evenly across the lobby. This uniform distribution makes the fractional representation an effective approximation. Expected values may therefore be used:

$$\mathbb{E}[L_K] \approx \frac{w_{V_F} \cdot \mathbb{E}[V_F] + w_{K_C} \cdot \mathbb{E}[K_C]}{w_{V_F} + w_{K_C}}$$

$D_P$ does not model accumulated knowledge explicitly, but it accounts for it implicitly. In theory, accumulation would follow a decaying rate. As more knowledge is acquired, the user has less left to learn, so the rate of learning new information falls, and $K_C$ approaches an asymptotic maximum. Because of this, the goal of keeping $D_P$ below a maximum value $D_{P,0}$ (Section 5.2) absorbs the behavior of accumulated knowledge. $D_{P,0}$ should be estimated after sufficiently long exposure, so that it already reflects the point where $K_C$ has largely saturated. Analyzing the inflated fraction directly is therefore a more concrete measurement of how detectable an inflated population is than modeling $K_C$ over time.

Note that the accumulation of knowledge $K_C$ is also countered by a complex factor regarding the decay of the value of knowledge. Notably, as time passes, knowledge would decay in value, with it eventually approaching a baseline minimum value based on unchangeable facts it provides. Any knowledge that provides information, on, say, the current state of the room in a given area at a given time would decay in value as time passed, and the effects of that information becomes largely ineffective to understanding the current state of the lobby. Understanding the rate of how each individual piece of knowledge decays is simply too complex, however, and this effect should also simply be absorbed into the measure of a saturated $K_C$. 

**Lobby chaos ($L_C$).** This is how visually busy the viewable part of the lobby is (or effectively how much information exists within a user's view), based on things like the number of objects and how much is moving. $L_C \in [0, 1]$, and its expected value $\mathbb{E}[L_C]$ is used. Chaos is the amount of information a user has to take in, so $L_C \uparrow \Rightarrow L_A \downarrow$. Notably, the player population $P$ acts as a source of chaos, particularly when analyzing whether a population is inflated. To recognize inflation, a user would have to count the real players, which becomes more difficult with more players. Chaos therefore grows with the number of players: $P \uparrow \Rightarrow \mathbb{E}[L_C] \uparrow$. 

Note that, we do not account for this effect in the final model. In general, we assert that this effect's weight on $D_P$ remains nearly constant. The reason for this is that, for a variety of reasons (performance, player density, mechanics, etc.) individual game lobbies almost never hold a large enough amount of players in order for the amount of players present to actually damage a user's ability to process information. Effectively, this model makes a reasonable player assumption, which does hold for the vast majority of games. If, however, there is an extremely large amount of players in a lobby (enough to hinder a user's ability to process information regarding these players), then this effect can still be absorbed in experimentation when estimating a $D_{P,0}$ threshold by properly simulating the large amount of players. Additionally, the effect would only decrease the true detectability of an inflated lobby, as the reasonable player model is the assumption made, which serves as a baseline. Therefore excluding this metric is conservative in nature, and avoids the complexity of modelling how the ability of a user to process information is hindered by more players producing more information, which is a metric that would require extra empirical analysis. Instead, absorbing this effect into $D_{P,0}$ serves as an effective way to handle the empirical results for this effect, while not actually requiring any extra empirical analysis. Note we will still use the relation of players to chaos to strengthen monotonicity logic, but the actual weight of this relation will be treated as negligible within the model. 

**User awareness ($U_A$).** This is how closely a user pays attention to detail, which decides how much of the available information the user can take in. $U_A \in [0, 1]$, and $U_A \uparrow \Rightarrow L_A \uparrow$. It is a hidden trait of each user that cannot be observed directly, so it is something of a nuisance variable. It depends on how willing and how able a user is to pay attention, so any estimate of it is likely to be both biased and inconsistent.

### 4.3 Reducing to two information effects

The three components affect two different properties of the information a user has. $L_K$ affects how *complete* that information is. $L_C$ and $U_A$ together affect how much of it can be *processed*: chaos is how much there is to process, and user awareness is how much a given user can handle. Other variables affect each of these as well, but the amount of available information is a strong factor shared by each. Because of this, the model can be reduced to two variables:

- $I_C \in [0, 1]$, information completeness
- $I_P \in [0, 1]$, information processability

Completeness goes up as lobby knowledge goes up:

$$\mathbb{E}[L_K] \uparrow \Rightarrow \mathbb{E}[I_C] \uparrow$$

Processability depends mainly on the information load $L_C$ and the user's capacity $U_A$:

$$\mathbb{E}[I_P] = \mathbb{E}\left[h(L_C, U_A)\right]$$

where $h$ is an unspecified function that decreases with $L_C$ and increases with $U_A$:

$$L_C \uparrow \Rightarrow I_P \downarrow, \qquad U_A \uparrow \Rightarrow I_P \uparrow$$

Since players are directly related to chaos, putting these together shows that processability drops as the lobby gains players, assuming $U_A$ remains constant:

$$P \uparrow \Rightarrow \mathbb{E}[L_C] \uparrow \Rightarrow \mathbb{E}[I_P] \downarrow$$

Both effects are inversely related to inflation. When the information is more complete or easier to process, inflated numbers are easier to spot, so less inflation can be used:

$$\mathbb{E}[I_C] \uparrow \Rightarrow F \downarrow, \qquad \mathbb{E}[I_P] \uparrow \Rightarrow F \downarrow$$

### 4.4 Total information awareness

To combine the two effects, total information awareness $I_A \in [0, 1]$ is defined as a weighted average, which serves as a measure of how well a user can piece information together:

$$I_A = \frac{w_{I_C} I_C + w_{I_P} I_P}{w_{I_C} + w_{I_P}}$$

where $w_{I_C}, w_{I_P} \geq 0$ and are not both zero. These weights set how much completeness and processability each factor in. Using expected values instead generalizes this analysis:

$$\mathbb{E}[I_A] = \frac{w_{I_C} \cdot \mathbb{E}[I_C] + w_{I_P} \cdot \mathbb{E}[I_P]}{w_{I_C} + w_{I_P}}$$

$I_A$ goes up with each of its inputs, or more specifically, $I_C \uparrow \Rightarrow I_A \uparrow$ and $I_P \uparrow \Rightarrow I_A \uparrow$, which follows from the partial derivatives:

$$\frac{\partial I_A}{\partial I_C} = \frac{w_{I_C}}{w_{I_C} + w_{I_P}} \geq 0, \qquad \frac{\partial I_A}{\partial I_P} = \frac{w_{I_P}}{w_{I_C} + w_{I_P}} \geq 0$$

Combining this with the relations above gives $L_C \uparrow \Rightarrow I_P \downarrow \Rightarrow I_A \downarrow$.

Recalling the relation of players to lobby chaos stated earlier:

$$P \uparrow \Rightarrow \mathbb{E}[L_C] \uparrow \Rightarrow \mathbb{E}[I_P] \downarrow \Rightarrow \mathbb{E}[I_A] \downarrow$$

### 4.5 What this means for $F$

A higher $\mathbb{E}[I_A]$ means users can piece information together better, which makes inflated numbers easier to detect, and hence $F$ would decrease. Therefore, $F$ depends on it inversely:

$$\mathbb{E}[I_A] \uparrow \Rightarrow F \downarrow, \qquad \text{equivalently} \qquad \mathbb{E}[I_A] \downarrow \Rightarrow F \uparrow$$

With $\mathbb{E}[L_K]$ held fixed, the two key relations are

$$P \uparrow \Rightarrow \mathbb{E}[I_A] \downarrow, \qquad \mathbb{E}[I_A] \downarrow \Rightarrow F \uparrow$$

By transitivity:

$$P \uparrow \Rightarrow F \uparrow$$

This matches the property found earlier (Section 3.2), so the information factors are directionally consistent with the proportional relationship between $F$ and $P$.

---

## 5. Form of $F(P)$

### 5.1 Growth constraint

Here $F$ is treated as a smooth function of $P$; whole-number rounding is introduced in Section 5.3. Since $D_P \approx \phi$, $\phi$ should not rise as $P$ grows, which means $\frac{F}{P}$ must not rise. If it did, detectability would increase, or more specifically, $\frac{F}{P} \uparrow \Rightarrow \phi \uparrow \Rightarrow D_P \uparrow$. Suppose $\frac{d^2 F}{dP^2} \gt 0$ (and knowing that $P > 0$ and that $F(0) = 0$, as applying inflation to an entirely empty lobby would be completely detectable). Then $\frac{F}{P}$ would rise with $P$, or $P \uparrow \Rightarrow \frac{F}{P} \uparrow$. We can prove this starting by differentiating the expression $\frac{F}{P}$ with respect to $P$, or hence taking $\frac{d}{dP}(\frac{F}{P})$:

<div align="center">

$$\frac{d}{dP}(\frac{F(P)}{P}) = \frac{P \cdot F'(P) - F(P)}{P^2}$$

</div>

We must show that $\frac{d}{dP}(\frac{F(P)}{P}) > 0$, meaning that $\frac{F}{P}$ would rise.

<div align="center">

$$\frac{P \cdot F'(P) - F(P)}{P^2} > 0 \Rightarrow P \cdot F'(P) - F(P) > 0$$

</div>

So we must show $P \cdot F'(P) - F(P) > 0$ in order to prove our claim that $P \uparrow \Rightarrow \frac{F}{P} \uparrow$.

To do this, we should differentiate $P \cdot F'(P) - F(P)$ with respect to $P$ and show the rate of change is greater than zero ($\frac{d}{dP}(P \cdot F'(P) - F(P)) > 0$), and then apply our initial condition of $F(0) = 0$ to show that the expression $P \cdot F'(P) - F(P) \ge 0$ for $P = 0$, meaning that it would always be positive if $\frac{d}{dP}(P \cdot F'(P) - F(P)) > 0$.

<div align="center">

$$\frac{d}{dP}(P \cdot F'(P) - F(P)) = P \cdot F''(P) + F'(P) - F'(P) > 0 \Rightarrow P \cdot F''(P) > 0$$

</div>

For $\frac{d^2 F}{dP^2} \gt 0$ and $P > 0$, this must be true.

This tells us that $P \cdot F'(P) - F(P)$ is increasing, and we know in order to prove our claim we must show that $P \cdot F'(P) - F(P) > 0$. Thus we need to show that for the minimum possible value of $P$, the expression $P \cdot F'(P) - F(P)$ is non-negative, or that for $P = 0$, then $P \cdot F'(P) - F(P) \ge 0$. If this is true, then for $P > 0$, then $P \uparrow \Rightarrow \frac{F}{P} \uparrow$.

To show $P \cdot F'(P) - F(P) \ge 0$ for $P = 0$, plug in $P = 0$ and $F(0) = 0$:

<div align="center">

$$0 \cdot F'(0) - F(0) = 0 - 0 = 0 \ge 0$$

</div>

therefore this is true. Thus, $P \uparrow \Rightarrow \frac{F}{P} \uparrow$.

Therefore $\phi$ will also increase, and therefore detectability would rise too, or more specifically, $P \uparrow \Rightarrow \frac{F}{P} \uparrow \Rightarrow \phi \uparrow \Rightarrow D_P \uparrow$. So $F$ can grow at most linearly in $P$, meaning $\frac{d^2 F}{dP^2} = 0$ (hence keeping detectability constant). Lower orders of growth, where $\frac{F}{P}$ shrinks, are still allowed.

### 5.2 The linear solution

Maximizing inflation at a given level of detectability calls for the fastest growth the constraint allows, which is linear. The inflated fraction does not change when both counts are scaled together. For any $c \gt 0$,

$$\phi(cP, cF) = \frac{cF}{cP + cF} = \frac{cF}{c(P+F)} = \frac{F}{P+F} = \phi(P, F)$$

so a linear $F$ keeps $D_P$ constant as $P$ grows. If $F(P) = kP$, then $F(cP) = cF(P)$, meaning $P$ and $F$ are scaled together by the same factor $c$, which was shown previously to keep $\phi$, and therefore $D_P$, constant.

Let $F(P) = kP$, where $k \gt 0$ is a free parameter. Then

$$D_P \Big|_{F=kP} = \frac{kP}{P + kP} = \frac{k}{k+1}$$

This does not depend on $P$. Note that $k$ is directly related to $D_P$, or $k \uparrow \Rightarrow D_P \uparrow$. Requiring $D_P \leq D_{P,0}$ therefore places a limit on $k$, such that $k \leq k_{\max}$. $k_{\max}$ is the value of $k$ where $D_P = D_{P,0}$, found by solving for $k$:

$$D_P = \frac{k}{k+1} \le D_{P,0}$$

$$\frac{k+1}{k} = 1 + \frac{1}{k} \ge \frac{1}{D_{P,0}}$$

$$\frac{1}{k} \ge \frac{1}{D_{P,0}} - 1$$

$$k \le \frac{1}{\frac{1}{D_{P,0}} - 1} = \frac{D_{P,0}}{1 - D_{P,0}} = k_{\max}$$

So $D_P \leq D_{P,0}$ holds exactly when $k$ is at or below this $k_{\max}$ cap (with $k \gt 0$ additionally enforced, since at $k = 0$ no inflation would occur):

$$0 \lt k \leq \frac{D_{P,0}}{1 - D_{P,0}}$$

### 5.3 Whole-number output

$F$ has to be a whole number, so rounding is needed. Though any rounding operation is allowable, the floor function is used here to minimize detectability:

$$F(P) = \lfloor kP \rfloor$$

Flooring can only lower $\phi$, so the bound on $D_P$ still holds, or more specifically, $\lfloor kP \rfloor \leq kP \Rightarrow \frac{\lfloor kP \rfloor}{P + \lfloor kP \rfloor} \leq \frac{k}{k+1}$. The system becomes:

$$F(P) = \lfloor kP \rfloor, \qquad D_P = \frac{k}{k+1} \leq D_{P,0}, \qquad 0 \lt k \leq \frac{D_{P,0}}{1 - D_{P,0}}$$

$D_{P,0}$ serves as an empirical value, and $k$ can be set within the provided constraint depending on $D_{P,0}$.

Note that this analysis covers the system at the steady state of inflation, meaning when no new inflation is being applied to the population. If the real population changes slightly, causing inflation to change slightly, that should still not act as a heavily disruptive transient, assuming a reasonable $\frac{dP}{dt}$. When inflation has just been newly applied, however, a system for time should be introduced, which is covered next.

---

## 6. Form of $F(P, t)$

### 6.1 Why time matters

If the full inflation is applied all at once, the displayed player count jumps suddenly. That jump is easy to notice. Slowing down how fast $F$ changes reduces this, or more specifically, $\left|\frac{\partial F}{\partial t}\right| \downarrow \Rightarrow D_T \downarrow$ (Section 6.2). Note that the absolute value of this metric should be measured: if inflation dropped (meaning $\frac{\partial F}{\partial t} < 0$), a user would also be more likely to recognize that inflation had been applied in the first place. This means removing inflation should use the same rate as applying inflation in order to keep the same temporal detectability. The elapsed time $t$ should therefore act as a damping factor that turns $F(P)$ into $F(P, t)$.

### 6.2 Temporal detectability $D_T$

$D_T$ should therefore depend on the magnitude of the rate of change of inflation, $\left|\frac{\partial F}{\partial t}\right|$. Additionally, the number of players $P$ is also a factor. For a given rate of change, the same change in $F$ is more significant in a smaller lobby than in a larger one, as larger populations (or lobbies) also tend to have larger fluctuations or variance, so greater rates of change in larger lobbies are expected, and vice versa for smaller lobbies. The fractional relation between $\left|\frac{\partial F}{\partial t}\right|$ and $P$ can therefore approximate $D_T$:

$$D_T \approx \frac{\left|\frac{\partial F}{\partial t}\right|}{P}$$

Similarly to population detectability, this is only an approximation. It also assumes that natural fluctuations scale roughly in proportion to $P$ (a constant relative volatility); the argument above establishes only that they grow with $P$. Other, more complex factors can drive $D_T$, the most notable and important of which is unusual influx.

### 6.3 Unusual influx

In reality, $D_T$ must hold some explicit relation to $dP/dt$, or more specifically, some fractional difference of $\frac{\mathbb{E}[dP/dt]}{dP/dt}$. Note here that $\mathbb{E}[dP/dt]$ would be expected to be an approximately periodic function, following the cycles by which players join and leave. If integrated over a single player cycle or period (for example, say a daily cycle of players who log on at approximately the same time every day), the total change of players would be expected to be small, or hence $\int_{T} \mathbb{E}[dP/dt] dt \approx 0$. This means the game is not gaining or losing any new players; it is only cycling the ones it already has. With daily player cycles, for example, a CCU may fluctuate over time daily, but the pattern of each daily cycle will generally remain approximately the same.

However, if external factors unusually drive players into or out of a game more quickly, then $D_T$ would in theory vary. For example, if a game received an update, or a famous social media creator made content on a game, it may drive an unusual number of players into the game quickly. Inflating more quickly during this period would make the inflation less noticeable.

It is important to note, however, that this effect applies the same to the reverse case, rather than oppositely. If, for some reason, many players left a game very quickly, inflating the player count would only decrease the rate of change of $P + F$ (the entire CCU). In theory, this would make players view the displayed rate of change as closer to $\mathbb{E}[dP/dt]$, because players would not necessarily know that the real player count was decreasing. Hence, in this special case where players drop, increasing the rate of change of the inflation would counteract the drop's effect on the total CCU, which would make players perceive a normal $dP/dt$ without knowledge of the drop, thus lowering detectability and user awareness. Effectively, this is a special effect of the reverse case of inflation, where players dropping too quickly could raise temporal detectability. This mirrors the idea that removing an inflation instantly could alarm players that an inflation occurred in the first place (Section 6.1). The difference is that, here, the uncontrollable variable ($P$) is the one rapidly dropping. Therefore, the rate of application of the controllable variable ($F$) can be increased to counteract this effect. However, it is critical to note that this concept does not apply to population detectability, as that is largely an effect of analyzing the actual game lobbies themselves. Because of this, it cannot be used to increase $F$ itself, but it can be used to increase the rate of its application over time, or $\left|\frac{\partial F}{\partial t}\right|$.

Notably, the masking effect captured by this ratio is absolute in manner, as it is the same in both directions, meaning:

$|dP/dt| \uparrow \Rightarrow \left|\frac{\mathbb{E}[dP/dt]}{dP/dt}\right| \downarrow \Rightarrow D_T \downarrow$.

This effect is rather complicated to predict and apply directly to a formula, however. To eliminate it, $dP/dt \approx \mathbb{E}[dP/dt]$ is assumed, removing the factor. This is especially acceptable as unusual influx would usually lower temporal detectability, leading to conservative measures. Regular player fluctuations would almost never decrease in magnitude, as the expected player fluctuations act as a base state, from which a disturbance causes excitement and thus an increase in the magnitude of the rate of change of players. This conservatism does not extend to a sudden mass departure of players: because $F$ is tied to $P$, a rapid drop in $P$ forces a rapid drop in $F$ and raises $D_T$, as described above. The model therefore assumes no such event. Handling one would mean trading a temporary rise in $D_P$ (keeping $F$ high while $P$ falls) against a rise in $D_T$ (letting $F$ fall with $P$).

### 6.4 The damping factor

Here we use a linear scaling for the damping factor. The time scaling factor rises linearly from $0$ to $1$ over $[0, T_0]$ and stays at $1$ afterward, where $t$ is the time since inflation began and $T_0$ is the damping duration:

$$
F(P, t) =
\begin{cases}
\left\lfloor kP \cdot \dfrac{t}{T_0} \right\rfloor & t \lt T_0 \\
\lfloor kP \rfloor & t \geq T_0
\end{cases}
$$

Equivalently, $F(P, t) = \left\lfloor kP \cdot \min\left(\frac{t}{T_0}, 1\right) \right\rfloor$.

### 6.5 Detectability during the damping period

Holding $P$ constant and differentiating without the floor operation (i.e., differentiating the continuous form so that the function is actually differentiable), for $t \lt T_0$:

$$\frac{\partial F}{\partial t} = \frac{kP}{T_0}$$

which does not change with $t$ and is positive, so $\left|\frac{\partial F}{\partial t}\right| = \frac{kP}{T_0}$. Substituting into $D_T$:

$$D_T = \frac{\frac{kP}{T_0}}{P} = \frac{k}{T_0} \qquad (t \lt T_0)$$

Note that $t$ and $T_0$ are measured in one standardized unit of time, the tick (for example, one minute), and are treated as plain numbers counted in ticks. $D_T$ is therefore a rate per tick, so a value of $D_{T,0}$ is only meaningful together with the tick length it was measured in, and the simulation's tick must match that length. For $t \geq T_0$ the time scaling factor is exactly $1$, so $\left|\frac{\partial F}{\partial t}\right| = 0$ and $D_T = 0$. The steady-state system analysis applies here (Section 5.3), so the temporal constraint only needs to be enforced during the damping period.

$T_0$ is assumed to be a discrete value, such as a whole number of minutes, so $T_0 \geq 1$. Additionally, $k$ is presumed not to exceed $1$. This is not a technical constraint of the system, but a general observation. $k = 1$ corresponds to $D_P = \frac{1}{2}$, i.e., doubling the real player count, which would almost certainly be detected. A larger $k$ would therefore never be chosen. Consequently, the greatest value of $D_T$ is $1/1 = 1$ and the lowest is $0$, reached once damping ends. Realistically, $D_T \in [0, 1]$ can therefore be bounded in this analysis (though the bound is not technically required, it will almost certainly hold).

The damping duration $T_0$ is inversely related to the rate of change of $F$ with time, and thus to $D_T$ by extension. More specifically, $T_0 \downarrow \Rightarrow \left|\frac{\partial F}{\partial t}\right| \uparrow \Rightarrow D_T \uparrow$. Because of this, requiring $D_T \leq D_{T,0}$ results in a minimum bound for $T_0$, such that $T_0 \ge T_{0,\min}$. $T_{0,\min}$ is the $T_0$ at which $D_T = D_{T,0}$, found by solving for $T_0$:

$$D_T = \frac{k}{T_0} \le D_{T,0}$$

$$\frac{T_0}{k} \ge \frac{1}{D_{T,0}}$$

$$T_0 \ge \frac{k}{D_{T,0}} = T_{0,\min}$$

Note that this constraint also requires knowledge of $k$, but since a constraint for $k$ was already defined, $k$ can be identified first, and then this constraint can be solved.

With these functions and inequalities, a complete population- and time-based system can now be pieced together.

---

## 7. Final System

### 7.1 Steady-state system

$$F(P) = \lfloor kP \rfloor, \qquad D_P = \frac{k}{k+1} \leq D_{P,0}, \qquad 0 \lt k \leq \frac{D_{P,0}}{1 - D_{P,0}}$$

**Empirical parameters:** $D_{P,0}$

**Constrained parameters:** $k$

### 7.2 Transient system (complete system)

$$
F(P, t) =
\begin{cases}
\left\lfloor kP \cdot \dfrac{t}{T_0} \right\rfloor & t \lt T_0 \\
\lfloor kP \rfloor & t \geq T_0
\end{cases}
$$

$$D_P = \frac{k}{k+1} \leq D_{P,0}$$

$$D_T = \frac{k}{T_0} \qquad (t \lt T_0)$$

$$D_T = 0 \qquad (t \geq T_0)$$

$$0 \lt k \leq \frac{D_{P,0}}{1 - D_{P,0}}, \qquad T_0 \ge \frac{k}{D_{T,0}} = T_{0,\min}$$

**Empirical parameters:** $D_{P,0}$, $D_{T,0}$

**Constrained parameters:** $k$, $T_0$

---

## 8. Extension to Multiple Lobbies with an Aggregate CCU

### 8.1 Superlobby definition

Now, suppose a game has $N \geq 2$ different lobbies. A CCU display is given not only for each individual lobby, but also for the sum of all lobbies. This sum is termed a "superlobby." This means

$$P_{total} + F_{total} = \sum_{n=1}^{N} (P_n + F_n)$$

where $P_{total} = \sum_{n=1}^{N} P_n$ and $F_{total} = \sum_{n=1}^{N} F_n$.

### 8.2 Nonexistence of $D_P$ for superlobbies

$D_P$ cannot exist for superlobbies. A superlobby, really, is just a CCU measure summing the CCUs of all individual lobbies together. Therefore, there is no actual game lobby in a superlobby. This means, by extension, there is also no $D_P$ as population detectability relies on the detection of inflation based on players acquiring information within the game lobby itself. Without an actual game lobby, there is no $D_P$ metric. Note, however, that this of course does not negate the individual $D_P$ values for each individual lobby, it only eliminates the existence of the $D_P$ value for a superlobby. 

### 8.3 Temporal detectability remains the same

It is important to note that the inflation analysis does not change for these superlobby systems, assuming the same $k$ and $T_0$ are applied to all servers, $D_T$ was found to be a constant value in the system. Therefore, it does not depend upon the current population. This means that whether inflation is applied to an individual lobby or a superlobby, the variation in population does not change the temporal detectability. As long as $k$ and $T_0$ remain fixed, $D_T$ is constant.

### 8.4 Smallest lobby requirement for $D_P$ analysis

As we previously found, there is no $D_P$ value for the superlobby itself, but, $D_P$ still exists for each individual lobby (Section 8.2). We can estimate a conservative $D_{P,0}$ for all individual lobbies by analyzing the single, most maximized $D_P$ lobby. Assuming the individual lobby counts are displayed, $D_{P,0}$ should be estimated using the smallest lobby. Recall that $\mathbb{E}[I_A]$ falls as $P$ grows, or that $P \uparrow \Rightarrow \mathbb{E}[I_A] \downarrow$ (Section 4.4). The consequence to this is that, while it is unaccounted for within the final system described, as the effect of it is simply too difficult to measure, larger lobbies would in theory pose a smaller population detectability, though this effect is likely negligible in most cases. However, because of this effect, it would in theory be best to analyze the smallest lobby to estimate the most extreme, or lowest value for $D_{P,0}$, which makes all other lobbies conservative in estimation and therefore the possibility of detection is negated. 

---

## 9. Simulations for Experimentation

### 9.1 Simulation justification

Now that a system describing inflation has been designed, the next requirement is to estimate empirical values. These purely empirical values are $D_{P,0}$ and $D_{T,0}$. Note that, while $k$ and $T_0$ are unknowns, they are bounded to a constraint based on the previously mentioned empirical values. In order to maximize inflation while staying below a certain detectability, $k$ and $T_0$ simply should be set to the most extreme condition possible along their bounds (that being $k_{max}$ and $T_{0,min}$ respectively, which prioritizes the final amount of inflation over how quickly it is reached, since $T_{0,min}$ grows with $k$), making these fully solvable based on the empirical values of $D_{P,0}$ and $D_{T,0}$. Hence, only $D_{P,0}$ and $D_{T,0}$ need to be determined experimentally. Simulating models similar to the real environments inflation may be applied to is a practical way to produce the controlled conditions an experiment needs to estimate proper $D_{P,0}$ and $D_{T,0}$ values. This paper does not design that experiment itself (its subjects, procedure, or analysis). Instead, it designs a simulation such an experiment could use, programmed in Rust, for the temporal detectability threshold case ($D_{T,0}$). The population detectability case ($D_{P,0}$) depends too heavily on the specific game to be simulated generically, as explained in Section 9.2. Effectively, theory for how to properly model and design a temporal simulation comparable to the environment of an actual game (at least to yield usable detectability values), while having far less complexity than an actual game, will be provided. 

### 9.2 Impracticality of modeling a population simulation

Building a simulation that can support accurate estimates of population detectability for a specific game is a much more complex task than doing so for temporal detectability, as precise modeling of the lobby itself is needed. When analyzing temporal detectability, the only non-negligible weight is the CCU counter itself, or more specifically how it changes with respect to time. With population detectability, however, how the lobby itself is composed must be modeled, so that the lobby knowledge a player acquires can be measured accurately in the simulation. 

Therefore, in order to simulate population detectability, an experimenter should effectively produce a lobby of their target game (assuming they have the source code) and simply get real players or advanced bots to play in these experimental lobbies. Then, they simply inflate the player count reading itself and produce a specific procedure to determine at what level the inflation is noticed by a subject. 

Essentially, modeling a simulation for population detectability is specific per game in order for the model's resolution to be fine-grained enough. Because of this, making a generic simulation is unnecessary for this. 

### 9.3 Why a temporal simulation can be modeled

Designing a temporal simulation, however, is much more justified. Analyzing temporal detection does not require physical lobby knowledge, meaning the entire world system required for the population simulation is unnecessary. Including the physical players of the lobby in this simulation's output would actually weaken any experiment built on it, since it would let subjects reach inflated-count answers by analyzing populations too. Instead, the only output required is the CCU over time, for enough ticks for a subject to sample. This makes the modeling for a temporal simulation almost general to any game, hence meaning designing a general simulation for temporal detectability is feasible. 

---

## 10. Temporal Simulation Design

### 10.1 Temporal simulation general rules

Since temporal detectability is unchanged for a superlobby (Section 8.3), we only require analysis of a singular equivalent lobby (whether the system being analyzed truly is a single lobby, or a superlobby). 

The population $P_{total}$ should start at the same constant value $P_0$ within each set of trials, mirroring the constant $P$ of the population simulation. Ideally, our $P_0$ should simply be estimated from actual CCU data of the target game as the average CCU of said game. Note, though, that $P_{total}$ may still fluctuate around $P_0$ to account for natural player variation, which is used only in this temporal simulation. $k$ should also be held constant within each set of trials, so that $T_0$ is the only inflation parameter that varies. The $k$ value should be taken from the one that is calculated to be the result of the population simulation. The $D_{P,0}$ estimated using that simulation is used to set $k$ to the largest value the constraint allows, which is $k = k_{\max} = \frac{D_{P,0}}{1 - D_{P,0}}$ (Section 5.2). This keeps the temporal trials consistent with the steady-state system that is actually intended for use.

Only the current CCU needs to be output for each tick, while a natural player variation function is still used to model natural CCU fluctuations. Inflation is then applied using a varying $T_0$. In trials where inflation does occur, the simulation should leave a sufficient period of time before or after the inflation starts, so that a subject in an experiment could effectively identify when the inflation actually occurred.

The simulation should also generate control trials, in which no inflation is applied at all, alongside the inflated trials. To keep subjects in an experiment from using position or order to guess which trials are controls, all trials should be randomly shuffled by the simulation before output.

### 10.2 Modeling natural player variation

A discrete player fluctuation model $P_f[t]$ can be defined such that $P[t + 1] = P[t] + P_f[t + 1]$.

Ideally, this is done by creating a probability distribution to select from, based on the specific data of a target game.

Here, the simulation's parameters must come from data collected from the actual game one would theoretically want to inflate or analyze the inflation potential of. CCU metrics should be collected and tracked over a sufficient amount of time to produce a large number of sample points based on a discretized time unit $t$. Minutes are probably the preferred choice of time unit. In general, the distribution of these fluctuations can be approximated as normal (which should be checked against the collected data), and its mean should be approximately 0. If the mean were nonzero to a non-negligible degree, the game would be actively growing or shrinking, which would violate the unusual influx assumption (Section 6.3) that the general change in players is equivalent to what would usually be expected for a steady player base that is neither growing nor shrinking. 

$\mathbb{E}[P]$ serves as the average and expected CCU. This becomes $P_0$ in the model.

For each point, the difference in players between it and the next point, $P[t + 1] - P[t]$, should be measured. From these differences, an average expected value along with its standard deviation should be computed to form a probability curve. Since $\mathbb{E}[P]$ is the stable average $P_0$, the mean $\mu$ of these differences should be approximately $0$ (a nonzero $\mu$ would make $P$ drift steadily away from $P_0$). $P_f[t]$ should randomly select a value per time tick based on this probability distribution. The result should be approximately normal, meaning

$P_f[t] \sim \mathcal{N}(\mu, \sigma^2)$, whose density is denoted $P_d(x)$.

In a program, a random selection that uses this curve as a probability weight can then be applied. The easiest way to do this is Monte Carlo rejection sampling (hit-and-miss), which proceeds by guessing and checking. A random point $x$ is selected within the probability curve, and $P_d(x)$ is found for that value. A random value $y \in [0, P_d(x)_{max}]$ is then generated (the maximum height of the normal curve is $\frac{1}{\sigma \sqrt{2\pi}}$), and the program checks whether $y \le P_d(x)$. If it is, the point is accepted and used as output. Otherwise, it is rejected and the process is retried until a point satisfies the conditions.

This technique can be used to select a random, unbiased $P_f[t]$ for each discrete time unit $t$. Note this makes $P_f[t]$ independent of the tick index $t$ (apart from the bounds introduced next), as fluctuations are instead selected randomly, weighted by a probability curve.

To prevent $P$ from straying far from $P_0$, one more addition can be made: minimum and maximum player bounds, which ensure $P$ stays within a range. The players will be bounded within the expected minimum and maximum daily players, which is also information that must be acquired from the target game and supplied to the simulation.

### 10.3 Temporal model specifics

The simulation outputs a basic set of trial vectors, which can be provided to subjects in an experiment. Each index represents a time unit, and each value within the indices represents the displayed aggregate CCU at that time ($P_{total} + F_{total}$).

The simulation should generate all trial vectors at the same size to reduce bias (if sizes varied, subjects in an experiment could begin to make heuristic assumptions about $T_0$ values). A general rule is that there must be at least as many time indices outside the transient inflation period as inside it, so that subjects could actively determine where the transient inflation starts and ends. Additionally, these non-transient sections can be divided between before the transient period (where no inflation has yet been applied) and after it (where the system is in steady-state inflation). The transient period can begin at any tick of the trial, from the very first tick up to the latest tick at which it still completes within the trial. Since ticks are indexed $0$ through $M - 1$ and the inflation first reaches its full value at tick $\text{start} + T_0$, that latest start is tick $M - T_0 - 1$. Because the simulation records exactly when each transient period starts and ends, if a subject in an experiment believes an inflation occurred, the subject can be asked for the specific period in which it occurred, and an error calculation can be made against the recorded ticks to analyze subject credibility. Control trials, which have no inflation applied, retain this same length to reduce bias, so that vector length never reveals whether a trial is a control.

Because all trial vectors should be the same size, the size $M$ must be determined by identifying the maximum number of ticks an inflation should take. In theory, this is the number of ticks the transient state consists of in the largest case, which is simply $T_{max}$, multiplied by 2 in order to enforce the general rule that the time outside the transient period should be at least that of the transient period:

$$M = 2 \cdot T_{max}$$

A largest tested value $T_{max}$ must therefore be selected for the trials. We can also set a smallest tested value $T_{min}$ to cut out values of $T_0$ that may be too low even by eye, and so do not need to be tested. These two bounds are chosen by the experimenter. They are distinct from the model's constraint $T_{0,\min} = \frac{k}{D_{T,0}}$ (Section 6.5), which is the value the experiment is meant to estimate, so $T_{min}$ should sit below it and $T_{max}$ above it. Along with this, a starting aggregate population $P_0$ and minimum and maximum player caps $P_{\min} \le P_0 \le P_{\max}$ must be selected for each set of trials, and $k$ is fixed to $k_{\max}$ as determined from the $D_{P,0}$ estimated using the population simulation. Here, $P_0$ represents the original total population.

The discrete natural player fluctuation model $P_f[t]$ is also applied per tick so that real player fluctuation is modeled, using the normal distribution (Section 10.2). Its parameters $\mu$ and $\sigma$ (and the minimum and maximum player caps) should be held constant within a trial set so that the fluctuation behavior does not become an additional varying factor. This leaves $D_T = \frac{k}{T_0}$ as a constant in each individual trial, but of course varies between trials.

The inflated trials go from $T_0 = T_{min}$ up to $T_0 = T_{max}$, and consecutive values are separated by a constant step $T_{jump}$. The $T_0$ values used are therefore $T_{min} + i \cdot T_{jump}$ for $i = 0, 1, 2, \dots$, for as long as $T_{min} + i \cdot T_{jump}$ does not exceed $T_{max}$, giving:

$$N_{inflated} = \left\lfloor \frac{T_{max} - T_{min}}{T_{jump}} \right\rfloor + 1$$

inflated trials. If $T_{max} - T_{min}$ is not a multiple of $T_{jump}$, the largest $T_0$ used falls short of $T_{max}$. Notably, this does not hurt the ratio of non-transient to transient ticks, since $M = 2 \cdot T_{max}$ is set from the upper bound of $T_{max}$ and therefore still leaves at least as many ticks outside the transient period as inside it (more in the scenario where $T_{0,largest} < T_{max}$). Control trials with no inflation at all are then added. The number of control trials to add is a design choice, set as a fraction $\rho_C$ of the number of inflated trials: $N_{control} = \text{round}(\rho_C \cdot N_{inflated})$. All of the trials are also scrambled randomly by the simulation, so that they can be handed to subjects in an experiment without their order revealing anything.

### 10.4 Programming the temporal simulation

The `rand = "0.8"` crate is used.

The program requires normal curve results from sampling, plus $P_0$, $k$ (which in actual use is not arbitrary, but is set to $k_{\max} = \frac{D_{P,0}}{1 - D_{P,0}}$ using the $D_{P,0}$ estimated with the population simulation), a $T_{min}$ and $T_{max}$ (`MIN_T` and `MAX_T`) to bound the tested range of $T_0$, a step $T_{jump}$ between consecutive $T_0$ values, and the control fraction. `MU` is set to `0.0` so that $P$ does not tend to drift away from $P_0$, and minimum and maximum player bounds are also implemented. `CONTROL_FRACT` is the number of control trials as a fraction of the number of inflated trials. `T_JUMP` is the increment between consecutive $T_0$ values across the inflated trials. It is a `usize` because it is passed directly to the iterator's `step_by`. The constants below are arbitrary dummy values, and are not necessarily the values that would be used in an experiment:

```rust
use rand::seq::SliceRandom;
use rand::Rng;
use std::collections::HashSet;
use std::f32;
use std::fs::File;
use std::io;
use std::io::{BufWriter, Write};

// Set constants.
const K: f32 = 0.2;
const MU: f32 = 0.0;
const SIG: f32 = 4.5; 
const P_0: i32 = 470;
const P_MIN: i32 = 350;
const P_MAX: i32 = 630;
const MIN_T: i32 = 1;
const MAX_T: i32 = 500;
const T_JUMP: usize = 1;
const CONTROL_FRACT: f32 = 1.0;
```

The $P_f[t]$ math should be implemented first. Assume some normal curve function (a mean and a standard deviation) is given, which the program can then parse into a normal curve function. A function can then be written that, given $x$ as a parameter (and using the constants `MU` and `SIG` for $\mu$ and $\sigma$), outputs the density $P_d(x)$ (named `norm` in the code).

The normal function is $f(x) = \frac{1}{\sigma \sqrt{2\pi}} e^{-\frac{1}{2}\left(\frac{x - \mu}{\sigma}\right)^2}$. Applying this:

```rust
// Generate Norm(x) for a given x.
fn norm(x: f32) -> f32 {
    let amplitude = 1.0 / (SIG * (2.0 * f32::consts::PI).sqrt());
    let exponent = -0.5 * ((x - MU) / SIG).powi(2);
    amplitude * exponent.exp()
}
```

To select a random $x$ value along the curve, the selection should be made within three standard deviations of the mean (as normal curves theoretically extend out to infinity, and three standard deviations cover 99.7% of the curve, as given by the empirical rule).

```rust
let x_bounds: [f32; 2] = [MU - 3.0 * SIG, MU + 3.0 * SIG];
```

The $y$ bounds extend from the minimum (0) to the maximum height ($\frac{1}{\sigma \sqrt{2\pi}}$):

```rust
let y_bounds: [f32; 2] = [0.0, 1.0 / (SIG * (2.0 * f32::consts::PI).sqrt())];
```

To perform a Monte Carlo guess, a random $x$ is selected along the $x$ bounds and rounded to the nearest whole player (fluctuations are whole numbers, and `norm` is evaluated at the rounded value), and a random $y$ along the $y$ bounds. The program then checks whether $y \le norm(x)$ and, additionally, whether the bounding constraint $P_{\min} \le p + x \le P_{\max}$ holds. If both hold, the point is accepted; otherwise, another guess is made. This calls for a loop.

```rust
// Monte Carlo hit & miss implementation.
fn monte_carlo(p: i32) -> i32 {
    let mut rng = rand::thread_rng();

    let x_bounds: [f32; 2] = [MU - 3.0 * SIG, MU + 3.0 * SIG];
    let y_bounds: [f32; 2] = [0.0, 1.0 / (SIG * (2.0 * f32::consts::PI).sqrt())];

    // Loop selections until a selected y <= norm(x).
    loop {
        let x_continuous = rng.gen_range(x_bounds[0]..=x_bounds[1]);
        let x_rounded = x_continuous.round();

        // Ensure the guess will not fluctuate players below or above bounds.
        if p + (x_rounded as i32) < P_MIN || p + (x_rounded as i32) > P_MAX {
            continue;
        }

        let y = rng.gen_range(y_bounds[0]..=y_bounds[1]);

        if y <= norm(x_rounded) {
            return x_rounded as i32;
        }
    }
}
```

These Monte Carlo iterations run per vector list point (one player fluctuation is generated per tick), so a function should be created to generate them.

```rust
// Generate a vector of player fluctuations per tick.
fn generate_player_fluctuations(len: usize) -> Vec<i32> {
    let mut result = Vec::new();
    let mut p = P_0;
    // Calculate each guess and track the p value with fluctuations to ensure it will stay within bounds;
    // the Monte Carlo function itself checks bounds to filter out guesses.
    for _ in 0..len {
        let monte_carlo_result = monte_carlo(p);
        result.push(monte_carlo_result);
        p += monte_carlo_result;
    }
    result
}
```

This completes player fluctuations. An actual individual trial vector must now be generated. A varying $T_0$ is used to inflate each trial, so it is taken as input. $T_0$ is an `i32` parameter, encapsulating it as a whole number representing minutes. The function also takes a boolean `inflate`. If `inflate` is `true`, the trial is generated as an inflated trial. If it is `false`, the entire $F$ calculation is skipped and $F$ is simply `0` at every tick, producing a control trial whose displayed count is just the real count. For an inflated trial, its $D_T = \frac{k}{T_0}$ is also computed. Alongside the vector and $D_T$, the function also returns the tick at which the inflation starts and the tick at which it ends, so the output is the tuple `(Vec<i32>, f32, usize, usize)`. The transient period covers ticks `start..end`, and the start tick is chosen uniformly at random from `0` through `len - t_0 - 1`, so that the inflation always reaches its full value on a tick that exists within the vector. For a control trial there is no transient, so the returned $D_T$, start, and end are `0.0`, `0`, and `0`.

```rust
// Build a single trial.
fn build_trial_vector(len: usize, t_0: i32, fluctuations: Vec<i32>, inflate: bool) -> (Vec<i32>, f32, usize, usize) {
    let mut rng = rand::thread_rng();

    let mut p = P_0;
    let mut result = Vec::new();
    // The latest valid start is len - t_0 - 1, so the inflation reaches its full value before the vector ends.
    let inflation_period_start = rng.gen_range(0..(len - (t_0 as usize)));
    let inflation_period_end = inflation_period_start + (t_0 as usize);

    // Fluctuate p, calculate f, push total CCU as p + f.
    for t in 0..len {
        p += fluctuations[t];

        let f = if inflate && t >= inflation_period_start {
            let time_fract = f32::min(1.0, ((t - inflation_period_start) as f32) / (t_0 as f32));

            f32::floor(K * (p as f32) * time_fract) as i32
        } else {
            0
        };

        result.push(p + f);
    }

    if inflate {
        (
            result,
            K / (t_0 as f32),
            inflation_period_start,
            inflation_period_end,
        )
    } else {
        (result, 0.0, 0, 0)
    }
}
```

Each trial needs a unique ID so that subject answers can later be matched back to the answer data without the order of the trials giving anything away. One random 4-digit hexadecimal ID is generated per trial, and repeats are not allowed. There are 65536 possible IDs, which is more than will be used per trial set.

```rust
// Generate random u16 value to then map to hexadecimal IDs.
fn generate_ids(count: usize) -> Vec<String> {
    let mut rng = rand::thread_rng();
    let mut seen: HashSet<u16> = HashSet::new();
    let mut ids = Vec::new();

    // Generate a random u16, check if the value has already been inserted.
    // If it is fresh, then insert it.
    while ids.len() < count {
        let value: u16 = rng.gen();
        if seen.insert(value) {
            ids.push(format!("{:04x}", value));
        }
    }

    ids
}
```

Multiple different trials are now generated. First, the inflated trials are generated from the $T_0$ bounds, `MIN_T` to `MAX_T`, advancing $T_0$ by `T_JUMP` after each trial (with `T_JUMP = 1`, every whole $T_0$ in the range gets a trial). Then, the control trials are generated. Their count is the number of inflated trials just built (which depends on `T_JUMP`) multiplied by `CONTROL_FRACT`, and each is built by calling the build function with `inflate` set to `false`. The controls are simply pushed onto the end of the list, since the whole list is randomly scrambled afterward, so their position at this stage does not matter.

Once all trials exist, each trial is assigned a unique ID. The result is split into two lists that share these IDs:

- The answer data, `(id, D_T, start, end)`.
- The vector data, `(id, vector)`, which is then completely randomly shuffled.

```rust
// Generate a full set of individual trials from MIN_T to MAX_T
// (in steps of T_JUMP), plus a set of control trials.
fn generate_trials(len: usize) -> (Vec<(String, f32, usize, usize)>, Vec<(String, Vec<i32>)>) {
    let mut trials: Vec<(Vec<i32>, f32, usize, usize)> = Vec::new();

    // Generate trials from MIN_T to MAX_T, advancing t_0 by T_JUMP each time.
    for t_0 in (MIN_T..=MAX_T).step_by(T_JUMP) {
        let fluctuations = generate_player_fluctuations(len);
        trials.push(build_trial_vector(len, t_0, fluctuations, true));
    }

    // The number of inflated trials built depends on T_JUMP, so count them directly.
    let inflated_count = trials.len() as f32;
    let control_count = (inflated_count * CONTROL_FRACT).round() as i32;

    // Build control trials.
    for _ in 0..control_count {
        let fluctuations = generate_player_fluctuations(len);
        trials.push(build_trial_vector(len, MIN_T, fluctuations, false));
    }

    let ids = generate_ids(trials.len());

    let mut data = Vec::new();
    let mut vectors = Vec::new();

    // Combine (zip) the IDs with the trial data, then iterate over each
    // and push the analysis data to one, and the trial data to the other,
    // with IDs passed to both.
    for ((trial, d_t, start, end), id) in trials.into_iter().zip(ids.into_iter()) {
        data.push((id.clone(), d_t, start, end));
        vectors.push((id, trial));
    }

    // Randomly order the trials to reduce bias.
    let mut rng = rand::thread_rng();
    vectors.shuffle(&mut rng);

    (data, vectors)
}
```

Each vector should now be output to a file, `./trials.txt`, with each vector on its own line. Each line begins with its position in the shuffled list, followed by the trial's ID in parentheses, then a period and the vector. For example:

```
1 (a3f9). [470, 473, 471, 468, ...]
2 (07bc). [470, 466, 469, 471, ...]
3 (e1d4). [470, 472, 475, 474, ...]
```

Each real line contains all $M$ values; they are truncated here.

These lines can then be sent to subjects in an experiment, with whatever instructions that experiment calls for. Only the shuffled vectors and their IDs are written to this file, so the file can be given to subjects without exposing the answer data. Because the order is random, the position number carries no information about $D_T$ or about whether a trial is a control.

```rust
// Write each trial in form:
// [number] (id). [...data]
// to ./trials.txt
fn write_trials(trials: &Vec<(String, Vec<i32>)>) -> io::Result<()> {
    let file = File::create("./trials.txt")?;
    let mut writer = BufWriter::new(file);

    for (i, (id, trial)) in trials.iter().enumerate() {
        let values: Vec<String> = trial.iter().map(|v| v.to_string()).collect();
        writeln!(writer, "{} ({}). [{}]", i + 1, id, values.join(", "))?;
    }

    writer.flush()?;
    Ok(())
}
```

The remaining values are written to a second file, `./trial_data.txt`. Each line uses the same ID format as `./trials.txt`, but the lines are written in generation order rather than shuffled order, so rows are matched to trials by ID rather than by position number. Each line consists of a line number, the ID in parentheses, then a period, followed by the $D_T$, the inflation start tick, and the inflation end tick, each separated by a space. For example:

```
1 (e1d4). 0.002 100 200
2 (a3f9). 0.0004 310 810
3 (07bc). 0 0 0
```

The last line is a control trial. The values shown are illustrative.

This file is kept separate from `./trials.txt` so that, in an experiment, it can be used to score subject answers (via the start and end ticks, and by checking which IDs are controls) and to compute the $D_T$ threshold (via $D_T$), without exposing any of it to the subjects.

```rust
// Write each trial's analysis data in form:
// [number] (id). [d_t] [start] [end]
// to ./trial_data.txt
fn write_trial_data(data: &Vec<(String, f32, usize, usize)>) -> io::Result<()> {
    let file = File::create("./trial_data.txt")?;
    let mut writer = BufWriter::new(file);

    for (i, (id, d_t, start, end)) in data.iter().enumerate() {
        writeln!(writer, "{} ({}). {} {} {}", i + 1, id, d_t, start, end)?;
    }

    writer.flush()?;
    Ok(())
}
```

Finally, everything can be put into a `main` function. Trials of size $M = 2 \cdot T_{max}$ (Section 10.3) are created, and then both `./trials.txt` and `./trial_data.txt` are written:

```rust
// Generate trials and write info to files.
fn main() -> io::Result<()> {
    let (data, vectors) = generate_trials((2 * MAX_T) as usize);
    write_trials(&vectors)?;
    write_trial_data(&data)
}
```

Altogether, this code can be pieced together into one complete simulation. See the temporal simulation files for the completed temporal simulation code.

---

## 11. Summary of Findings

### 11.1 System derivation

We have derived a system, for both a steady and transient state, for describing how to inflate a CCU by a number of fake players $F$, which in the steady state model acts as $F(P)$, and in the transient model acts as $F(P, t)$. We found that a linear model, in both the scaling factor of inflated to real players, along with the temporal damping factor, serve as effective fits to maximize the rate and amount of inflated players $F$, while still maintaining constant levels of detectability $D_P$ and $D_T$. We also extended this logic to "superlobby" systems, and found the resulting analysis would remain largely unchanged due to the nonexistence of $D_P$ for the superlobby itself, and the fact $D_T$ remains constant for the superlobby as it does not rely on the population $P$.

### 11.2 Empirical analysis

We determined that the variables $k$ and $T_0$ are constraint bound values, that, for the most effective results, could simply be set to their most extreme allowed values ($k_{\max}$ and $T_{0,\min}$), meaning their values are determinable and not empirical. Instead, the empirical values come in the form of $D_{P,0}$ and $D_{T,0}$, which serve as detectability threshold levels, and must be found via experimentation. 

### 11.3 Simulation design

We found that designing a simulation to describe $D_{P,0}$ would be infeasible, as a general simulation cannot hold the amount of resolution that specific individual games would require for this analysis. However, we determined this was not the case for developing a simulation to find $D_{T,0}$, and we therefore designed an example simulation that could be used in an experiment to do such.
