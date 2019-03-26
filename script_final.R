library(lme4)
library(lmerTest)
library(fitdistrplus)
library(emmeans)
library(tidyverse)

# read in the data
data <- read_csv("Mariesdata2.csv")

# tidy and wrangel the data to get it in the correct format
tidy_data <- data %>% 
          select(Participant, Item, type, 
          Focus, Cleft, YNResp.corr_raw, YNResp.rt_raw, changeResponse.corr_raw) %>%
          filter(type != "filler") %>% 
          select(-type)

colnames(tidy_data) <- c("Participant", "Item", "Focus", "Cleft", 
                         "Response", "RT", "Correct_word")

tidy_data[tidy_data$Focus == "y",]$Focus <- "Early"
tidy_data[tidy_data$Focus == "n",]$Focus <- "Late"

tidy_data$Focus <- as.factor(tidy_data$Focus)
tidy_data$Cleft <- as.factor(tidy_data$Cleft)
tidy_data$Participant <- as.factor(tidy_data$Participant)
tidy_data$Item <- as.factor(tidy_data$Item)

# set up contrasts for the two factors
my_contrasts <- matrix(c(.5, -.5))

contrasts(tidy_data$Focus) <- my_contrasts
contrasts(tidy_data$Cleft) <- my_contrasts

# convert RT variable to milliseconds
tidy_data$RT <- tidy_data$RT * 1000

# generate summary data and visualise
# all responses
tidy_data %>% 
  group_by(Focus, Cleft) %>%
  summarise(mean(Response), sd(Response), mean(RT, na.rm = TRUE), sd(RT, na.rm = TRUE))

# only when response is correct
tidy_data %>% 
  filter(Response == "1") %>%
  group_by(Focus, Cleft) %>%
  summarise(mean(Response), sd(Response), mean(RT, na.rm = TRUE), sd(RT, na.rm = TRUE))

# RT data when response is correct
tidy_data %>% filter(Response == "1") %>%
  ggplot(aes(x = Focus:Cleft, y = RT, colour = Focus:Cleft)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) 

# build binomial model on response data
model <- glmer(Response ~ Focus * Cleft + (1 | Participant) +
  (1 | Item), data = tidy_data, family = "binomial")
summary(model)

# binomial model on typed response data when people get it correct
corr_data <- tidy_data %>% filter(Response == "1")
model <- glmer(Correct_word ~ Focus * Cleft + (1 | Participant) +
                 (1 | Item), data = corr_data, family = "binomial")
summary(model)
emmeans(model, pairwise ~ Focus)
emmeans(model, pairwise ~ Cleft)

corr_data %>% 
  filter(Response == "1") %>%
  group_by(Focus, Cleft) %>%
  summarise(mean(Correct_word))

corr_data %>% filter(Response == "1") %>%
  ggplot(aes(x = Focus:Cleft, y = Correct_word, colour = Focus:Cleft)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) 

# main effect of focus - let's look at the descriptives
tidy_data %>% 
  group_by(Focus) %>%
  summarise(mean(Response), mean(RT, na.rm = TRUE))

# RT analysis only when people get the response right - first with gaussian distribution
# to find a model that converges and isn't over-parameterised (as indicated by
# the singular fit error) we drop the interaction term from the item random
# effect
correct_only_data <- filter(tidy_data, Response == "1")
model <- lmer(RT ~ Focus * Cleft + (1 + Focus + Cleft | Participant) +
                (1 + Focus * Cleft | Item), data = correct_only_data)

summary(model)
qqline(residuals(model))

# distribution of residuals above doesn't look normal so now try fitting 
# assuming sampling from the Gamma distribution.
# The data look like they follow the Gamma distribution:
descdist(correct_only_data$RT)
model <- glmer(RT ~ Focus * Cleft + (1 | Participant) +
                (1 | Item), data = correct_only_data, family = "Gamma")

# Even the intercept only model above does not converge
# try log transforming and then building linear model
correct_only_data <- correct_only_data %>% 
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Focus * Cleft + (1 + Focus + Cleft | Participant) +
                (1 + Focus + Cleft | Item), data = correct_only_data)
summary(model)
qqnorm(residuals(model))
qqline(residuals(model))

emmeans(model, pairwise ~ Focus * Cleft)

# RT analysis regarldess of whether people get the response correct
all_data <- tidy_data %>% 
  mutate(logRT = log(RT))
model <- lmer(logRT ~ Focus * Cleft + (1 + Focus | Participant) +
                (1 + Focus | Item), data = all_data)

summary(model)
qqnorm(residuals(model))
qqline(residuals(model))

emmeans(model, pairwise ~ Focus * Cleft)
  
# Plots speed vs accuracy
 tidy_data$Response <- as.factor((tidy_data$Response))

tidy_data %>% 
  ggplot(aes(x = Focus:Cleft, y = RT, colour = Response)) +
  geom_jitter(width = .2, alpha = .8, size = .5) +
  guides(colour = FALSE) + 
  scale_colour_hue(h = c(0, 360))

tidy_data %>% 
  ggplot(aes(x = Focus:Cleft:Response, y = RT, colour = Response)) +
  geom_jitter(width = .2, alpha = .8, size = .5) +
  guides(colour = FALSE) + 
  scale_colour_hue(h = c(0, 360)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black")

