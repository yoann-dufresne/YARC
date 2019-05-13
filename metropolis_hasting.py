from random import Random


class MetropolisHasting1D:
    """ Class implementing the metropolis hastings sub-sampling to mimic experimental distribution
    passed in input. To have a f function over all the distribution (even where there are no values)
    we connect successive (x, y) values by linear segments.

    For more reproducibility, the random seed used can be fixed when an object is created.
    """

    def __init__(self, distribution, num_steps, std_dev=None, seed=None, burning_steps=10000):
        """ At the begining the random uniform function can be initialize using the inputed seed.

        The other goal of the init function is to perform a default burn of the first values that
        can be outputed by the metropolis hastings algorithm. The number of burning step can be
        changed using the value of the parameter burning_steps.

        Arguments:
            distribution: An array containing all the tuples (x, y) as f(x)->y. This distribution
            is used to infer all the missing values between all pair of successive x. WARNING: the
            values must be sorted by increasing value of x.
            Here an example of distribution that can be passed:
            [(3,1), (4,1), (6,3), (7,5), (12, 2), (17,1)]

            std_dev: The sigma value for the gauss function from std python random lib. This will
            be used to find the next x value.

            num_step: Replace the step size argument by spliting the distribution into num_step
            different x values. The step size will be

            seed: An int value for initializing the random seed of the uniform random function. By
            default, the seed is choosen by the random class.

            burning_steps: An int value to indicate the number of step that has to be done before
            to consider the metropolis hastings as initilized.

        Exceptions:
            Exeception("Unordered distribution"): If the distribution is not ordered on the x value
        """

        # Verify that the distribution is correctly ordered and raise an exception if not
        previous = distribution[0][0]
        for x, _ in distribution[1:]:
            if x <= previous:
                raise Exception("Unordered distribution")
                return
            previous = x

        # fix std deviation if not set by user
        if not std_dev:
            std_dev = (distribution[-1][0] - distribution[0][0]) / 100.0
        self.std_dev = std_dev

        # # change all values to floats
        # distribution = [(float(x), float(y)) for x,y in distribution]

        step = (distribution[-1][0] - distribution[0][0]) / num_steps

        # Creating the distribution with the missing values
        self.distribution = []

        left_x_idx = -1
        right_x_idx = 0
        for x in self.__float_range(distribution[0][0], distribution[-1][0], step):
            while x > distribution[right_x_idx][0]:
                left_x_idx += 1
                right_x_idx += 1

            # If the value y exist, copy it
            if distribution[right_x_idx][0] == x:
                self.distribution.append(distribution[right_x_idx])
            # If y does not exist, do a linear interpolation
            else:
                self.distribution.append((x, distribution[left_x_idx][1] + \
                    (distribution[right_x_idx][1] - distribution[left_x_idx][1]) * \
                    (x - distribution[left_x_idx][0]) / \
                    (distribution[right_x_idx][0] - distribution[left_x_idx][0])
                ))

        # Initialize a random generator for this instance
        self.random = Random()
        if seed:
            self.random.seed(seed)

        # Compute the avg y value
        nb_values = 0
        total_value = 0.0
        for x, y in distribution:
            nb_values += y
            total_value += y * x
        
        avg_value = total_value / nb_values

        # Select closest x for average
        dy = abs(distribution[0][1] - avg_value)
        self.current_x, self.current_y = distribution[0]
        
        for x, y in self.distribution:
            cur_dy = abs(y - avg_value)

            if cur_dy < dy:
                dy = cur_dy
                self.current_x = x
                self.current_y = y

        # Burning step
        self.nb_accepted = 0
        self.nb_steps_performed = 0

        for _ in range(burning_steps):
            self.next_value()

        # Save burning acceptance
        self.burning_nb_accepted = self.nb_accepted
        self.burning_nb_steps_performed = self.nb_steps_performed
        self.nb_accepted = 0
        self.nb_steps_performed = 0


    def __float_range(self, start, stop, step):
        current = start

        while current <= stop:
            yield current
            current += step


    def next_value(self):
        """ Return the next value of the distribution using the metropolis hastings algorithm.
        The step start at the value self.current_x and go left or right regarding the standard
        deviation. Under min_x and over max_x, f(x) = 0

        Return float
            The float returned correspond to step value.
        """
        
        self.nb_steps_performed += 1

        # Get next values
        next_x = self.random.gauss(self.current_x, self.std_dev)
        next_y = 0
        # Compute next_y
        if next_x >= self.distribution[0][0] and next_x <= self.distribution[-1][0]:
            idx = (len(self.distribution)-1) * \
                (next_x - self.distribution[0][0]) / \
                (self.distribution[-1][0] - self.distribution[0][0])
            
            next_y = self.distribution[round(idx)][1]

        # Compute proba
        proba = min(1.0, next_y / self.current_y)

        if self.random.uniform(0, 1) < proba:
            self.current_x = next_x
            self.current_y = next_y

        return (self.current_x, self.current_y)


    def next_values(self, nb_val):
        return [self.next_value() for _ in range(nb_val)]


if __name__ == "__main__":
    distrib = [(3,1), (4,1), (6,3), (7,5), (12, 2), (17,1)]
    metro = MetropolisHasting1D(distrib, 2, burning_steps=10)

    for _ in range(100):
        print(metro.next_value())
