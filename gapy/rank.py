def rank_dominate_by_feature(individuals: list, weighting: list = None):

    """Ranks individuals by the number of individual features they dominate. The individual that dominates 
    the most other individual features corresponds to rank 1.

    Arguments:
        individuals (list[Individual]): List of individuals.
        weighting (list[float]): List of weights to specify the importance of features.

    Returns:
        list[int]: The list of rankings corresponding to the individuals.
    """

    domination_numbers = [0 for _ in individuals]
    for i, individual_1 in enumerate(individuals):
        for j, individual_2 in enumerate(individuals):

            if i == j:
                continue
            
            dominating_features = individual_1.get_dominating_features(individual_2)
            for k, feature in enumerate(dominating_features):
                if feature:
                    # consider weighting of features
                    if weighting is not None:
                        domination_numbers[i] += 1 * weighting[k]
                    else:
                        domination_numbers[i] += 1

    unique_domination_numbers = sorted(list(set(domination_numbers)), reverse=True)

    rankings = [None for _ in individuals]
    for i, unique_domination_number in enumerate(unique_domination_numbers):
        for j, _ in enumerate(domination_numbers):
            if domination_numbers[j] == unique_domination_number:
                rankings[j] = i + 1

    return rankings

def rank_is_dominated_by_feature(individuals: list, weighting: list = None):

    """Ranks individuals by the number of individual features they are dominated by. The individual that is dominated 
    by the least other individual features corresponds to rank 1.

    Arguments:
        individuals (list[Individual]): List of individuals.
        weighting (list[float]): List of weights to specify the importance of features.

    Returns:
        list[int]: The list of rankings corresponding to the individuals.
    """

    domination_numbers = [0 for _ in individuals]
    for i, individual_1 in enumerate(individuals):
        for j, individual_2 in enumerate(individuals):

            if i == j:
                continue
            
            dominating_features = individual_1.get_dominated_features(individual_2)
            for k, feature in enumerate(dominating_features):
                if feature:
                    # consider weighting of features
                    if weighting is not None:
                        domination_numbers[i] += 1 * weighting[k]
                    else:
                        domination_numbers[i] += 1

    unique_domination_numbers = sorted(list(set(domination_numbers)), reverse=False)

    rankings = [None for _ in individuals]
    for i, unique_domination_number in enumerate(unique_domination_numbers):
        for j, _ in enumerate(domination_numbers):
            if domination_numbers[j] == unique_domination_number:
                rankings[j] = i + 1

    return rankings


def rank_dominate(individuals: list):

    """Ranks individuals by the number of individuals they dominate. The individual that dominates 
    the most other individuals corresponds to rank 1.

    Arguments:
        individuals (list[Individual]): List of individuals.

    Returns:
        list[int]: The list of rankings corresponding to the individuals.
    """

    domination_numbers = [0 for _ in individuals]
    for i, individual_1 in enumerate(individuals):
        for j, individual_2 in enumerate(individuals):

            if i == j:
                continue
                
            if individual_1.dominates(individual_2):
                domination_numbers[i] += 1

    unique_domination_numbers = sorted(list(set(domination_numbers)), reverse=True)

    rankings = [None for _ in individuals]
    for i, unique_domination_number in enumerate(unique_domination_numbers):
        for j, _ in enumerate(domination_numbers):
            if domination_numbers[j] == unique_domination_number:
                rankings[j] = i + 1

    return rankings

def rank_is_dominated(individuals: list):

    """Ranks individuals by the number of individuals they are dominated by. The individual that is dominated
    by the least other individuals corresponds to rank 1.

    Arguments:
        individuals (list[Individual]): List of individuals.

    Returns:
        list[int]: The list of rankings corresponding to the individuals.
    """

    domination_numbers = [0 for _ in individuals]
    for i, individual_1 in enumerate(individuals):
        for j, individual_2 in enumerate(individuals):

            if i == j:
                continue
                
            if individual_1.is_dominated_by(individual_2):
                domination_numbers[i] += 1

    unique_domination_numbers = sorted(list(set(domination_numbers)), reverse=False)

    rankings = [None for _ in individuals]
    for i, unique_domination_number in enumerate(unique_domination_numbers):
        for j, _ in enumerate(domination_numbers):
            if domination_numbers[j] == unique_domination_number:
                rankings[j] = i + 1

    return rankings

def rank_non_dominated_fronts(individuals: list):

    """Ranks individuals by their membership to non dominated fronts.

    Arguments:
        individuals (list[Individual]): List of individuals.

    Returns:
        list[int]: The list of rankings corresponding to the individuals.
    """

    non_dominated_fronts = []
    individuals_idx = [i for i in range(len(individuals))]
    while len(individuals_idx) > 0:
        
        non_dominated_front = []
        for i in individuals_idx:
            domination_number = 0
            for j in individuals_idx:

                if i == j:
                    continue

                if individuals[i].is_dominated_by(individuals[j]):
                    domination_number += 1

            if domination_number == 0:
                non_dominated_front.append(i)

        non_dominated_fronts.append(non_dominated_front)

        for i in non_dominated_front:
            del individuals_idx[individuals_idx.index(i)]

    # assign ranks according to membership to non-dominated fronts

    rankings = [None for _ in individuals]
    for i, non_dominated_front in enumerate(non_dominated_fronts):
        for idx in non_dominated_front:
            rankings[idx] = i + 1

    return rankings