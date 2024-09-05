import math
import scipy.stats as stats
from statsmodels.stats.power import TTestIndPower

# Given data
n = 108  # total sample size
x = 48   # number of successes (sensors hit)
p = x / n  # observed proportion

# Calculate confidence interval
confidence_level = 0.95
z = stats.norm.ppf((1 + confidence_level) / 2)
margin_of_error = z * math.sqrt((p * (1 - p)) / n)

lower_bound = max(0, p - margin_of_error)
upper_bound = min(1, p + margin_of_error)

print(f"Observed proportion: {p:.4f}")
print(f"{confidence_level*100}% Confidence Interval: ({lower_bound:.4f}, {upper_bound:.4f})")

n = 108  # total sample size
observed_prop = 48 / 108  # observed proportion

# Define the effect size you want to detect
# For example, let's say you want to detect a 10% difference
null_prop = 0.5  # null hypothesis proportion
effect_size = abs(observed_prop - null_prop)

# Calculate the power
alpha = 0.05  # significance level
power_analysis = TTestIndPower()
power = power_analysis.solve_power(effect_size=effect_size, nobs1=n, alpha=alpha, ratio=1.0, alternative='two-sided')

print(f"Power: {power:.4f}")

# If you want to know the required sample size for a certain power
desired_power = 0.8
required_n = power_analysis.solve_power(effect_size=effect_size, power=desired_power, alpha=alpha, ratio=1.0, alternative='two-sided')

print(f"Required sample size for {desired_power*100}% power: {math.ceil(required_n)}")
