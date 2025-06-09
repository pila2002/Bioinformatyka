<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" class="logo" width="120"/>

# Elective Course: Bioinformatics Laboratory

**File Version:** v1.5, March 27, 2025
**Author:** Dr. Eng. Marcin Radom

---

## Project – Completion – General Rules

- Work in teams of two.
- **Main Report:**
    - Start with a description of the algorithm. Focus primarily on explaining the algorithm, whether it is my version from the provided materials, a modified one (explain what and where was changed), or an entirely different idea.
        - You can skip describing the instance generator unless you have made an interesting modification worth mentioning (which I doubt…).
        - Similarly, skip describing the graph creation module unless you made an original contribution in this part.
    - **Metaheuristics:** Provide a detailed description. For example, if using a genetic algorithm, explain exactly how crossover, mutation, and selection work. If it’s Tabu Search, describe how the tabu list works, whether there is an aspiration criterion, and how the algorithm deals with local optima, etc. Do the same for other algorithms, both standard and custom, so I can understand how the algorithm works without reading the code line by line.
    - **Tests and Results:** Present test results in tables and charts (see below for details; the quality of this section will largely determine your initial project assessment).
- Reports in PDF/DOC/ODF format can be sent by email, along with the complete project code (be careful not to include executable files in attachments). If I do not confirm receipt within 24-48 hours, assume the email was lost (usually due to spam filters or attachments).
- **Preferred Programming Languages:** C++, C\#, Java, Python.

---

## Project Completion Dates and Deadlines

- On Sunday, 15.06, groups L2 and L4 have priority, as these are their classes. I will be at work from about 9:00 to 16:00 or until the last scheduled group. If all L2 and L4 groups submit their reports on time (see below), discussions with them will take less than half the allotted time. Therefore, L1, L3, and L5 groups can come that day if they find time between other classes. Prior notification is appreciated (not required for L2, L4; they should only notify by email if they will not attend and want a Zoom meeting—see below).
- To ensure a smooth project discussion, I would like to receive the report by 9/10.06, i.e., by midnight or Tuesday morning, 10.06.
    - I will first reply by email to those who send their projects by this time to confirm receipt (if you don't get a reply, your email didn't arrive, as I reply to these as a priority). A day or two later, I will send another email with your scheduled time on Sunday to avoid unnecessary waiting in line.
    - If I receive your report between 10.06 and 12.06, it depends. If I manage to review it before Friday, you can come on Sunday.
- **Saturday, 14.06 – classes are canceled.**
<gothic mode>By decree of the benevolent Dean's Office, grades must be entered into USOS by 22.06.</gothic>
Therefore, we can have the "discussion" on Zoom. In this case, I must receive the report by Saturday, 14.06, at midnight at the latest. This way, on Sunday, instead of talking to you, I won't be bored at PP, but will be printing and reading. On Monday, 16.06, I will send emails with Zoom links. Honestly, I would prefer to unify this: for example, let's agree to hold Zoom discussions on Wednesday, 18.06, or Thursday, 19.06, around 18:00 or later. I simply don't want to have to schedule 15 Zoom meetings scattered across all days and hours.

---

## Testing

### Types of Tests

- **Metaheuristic Tuning Tests:** Tests for tuning metaheuristics.
- **DNA Sequencing Quality Tests:** Tests for the quality of DNA sequencing.


### Algorithm Parameter Testing

Metaheuristic (and advanced) algorithms contain numerous parameters that require specific values. For example:

- Genetic algorithm: mutation and crossover probabilities, selection rules, etc.
- Tabu Search: tabu list length, long-term memory, local optimum escape parameters.
- Ant Colony Optimization: number of ants, pheromone matrix management functions.

Before demonstrating how well your heuristics sequence DNA with various charts, justify these parameters with at least one or two simple charts, so they are less "dreamed up" and more justified. For example, if the genetic algorithm population size is 42, include a chart showing why 42 is optimal compared to 41 or 43. Also consider runtime, mutation/crossover probabilities, ant decisions, tabu list length, etc.

Do not overdo it—testing everything is impossible, even mathematically (multidimensional, multicriteria problem). But even one simple chart testing something from the algorithm will work wonders for my future complaints, or lack thereof, during review.

### Sequencing Tests

The main thing to test is how well the algorithm assembles DNA. There are three problems to choose from:

1. Classic SBH problem with negative and positive errors.
2. SBH problem with positional information and all types of errors.
3. SBH problem with repeat information and all types of errors.

**Testable Problem Parameters:**

- **DNA length (n):** Between 300 and 1000 (for classic SBH, 1000 is very high, especially with spectrum errors).
- **Oligonucleotide probe length (k):** Between 7 and 10. Larger values can be used but are less realistic—longer probes make assembly easier, shorter k increases negative errors from repeats.
- **Error types:**
    - Negative (from repeats—standard, assumed for longer DNA).
    - Negative "normal".
    - Positive.
- **Number/percentage of errors:**
2-3% of (n–k+1) is the absolute minimum. 10% negative errors is a lot, 15% is extreme. For positive errors, even 10-15% is a lot, but the algorithm may still handle it.
For combined cases (e.g., 10% negative and 10% positive), this is a high error rate and expect weaker results.
- **[Problem 2]** Test the precision of repeat information: typical info is "1 or more" or "1, 2, or more repeats". Exact counts (e.g., 7 repeats) are rare and would make things easier, but in practice, the two mentioned schemes are most common.
- **[Problem 3]** Parameter: range/precision of localization. Knowing an oligonucleotide is between positions 10 and 20 is much more helpful than between 10 and 150. Test how range size affects solution quality.

You do not need to test everything, but the more tests for different parameters, the better the report and final grade.

---

## Data Representation

Tables are fine as long as they are well-constructed and clearly labeled. Excel tables can be converted to charts. Key points for charts (based on Janina Bąk’s book “Statistically Speaking…”):

1. A chart is not a Valentine’s card—it must be labeled (title and axes).
2. All chart elements (colors, point shapes, etc.) must be appropriately named, so there is no doubt what they represent.
3. 25% of a large pizza is not the same as 25% of a small pizza—always include the total number of observations (N) covered by the chart.
4. Visualization is not exercise—make it easy for the audience to understand, e.g., by adding value labels.
5. A chart is not an Easter egg—use colors sparingly.
6. Be considerate of people who are overwhelmed by numbers—show only what is necessary and highlight what’s most important. Don’t overdecorate.

---

## Testing – In Other Words

The above has sufficed so far (more or less). The following is a more detailed repetition of the same testing idea, developed for a similar course in 2023-2024.

### Project Summary (as explained in the lab):

1. **Teamwork:** Work in pairs to implement algorithm(s) solving the sequencing by hybridization problem, including:
    - **Instance generator:** Generates hybridization spectrum, potentially with positive and negative errors (amount depends on the test). Typical DNA size: n > 300nt and n < 700nt, oligonucleotide size: 7 ≤ k ≤ 10.
    - **Initial solution generator:** Called the “random solution generator” in labs, as it reflects the usual “quality” of its output.
    - **Metaheuristic algorithm:** To score above 4.0–4.5, attempt to implement a metaheuristic that improves DNA reconstruction from the spectrum.
    - **Report:** Should include:
        - Description of created algorithms.
        - Algorithm parameter tests (choose 2–3 key ones for metaheuristics and test their impact on solution quality).
        - Problem instance tests (e.g., effect of hybridization error percentage on solution quality, effect of DNA/k length, etc.).
2. **For charts/tables,** use a string distance metric, e.g., Levenshtein distance. Since any problem instance starts with generated DNA, after reconstructing DNA from the (possibly erroneous) spectrum, compare the result to the original DNA to objectively assess the difference—these are your data points for charts.
3. **Classic SBH:** Known DNA length, oligonucleotide length (k), number of negative/positive errors, starting vertex.
4. **Test types in the report:**
    - **Type 1:** Tests the impact of algorithm parameters on solution quality.
    - **Type 2:** Tests the impact of problem instance variables (e.g., number of errors) on the difficulty of finding a good solution.

**Algorithm Parameter Testing Procedure:**
    - Use a fixed, representative set of test instances (at least several dozen, similar in k and n, differing in error types/amounts).
    - Example: To test the number of ants in ACO, run the metaheuristic for each value (e.g., 50, 100, 150, 200, 250 ants) on all instances, calculate the mean Levenshtein distance for each setting, and plot these as points on a chart.
    - Repeat for 2–3 parameters, even if there are more in the algorithm.
    - If you don’t use a metaheuristic, find any significant parameters in your naive algorithm and test them similarly.

**Problem Instance Parameter Testing Procedure:**
    - Generate a set of test instances (e.g., 20 DNA chains with identical n and k).
    - For each error percentage (e.g., 2%, 4%, 6%, 8%, 10%), modify the spectrum accordingly, run the algorithm, compare reconstructed DNA to the original, and calculate the mean Levenshtein distance.
    - Always use the same DNA chains for each error percentage to ensure only the error rate is being tested.

---

*End of translation.*

<div style="text-align: center">⁂</div>

[^1]: paste.txt

