library(lme4)
library(lmerTest)
library(fitdistrplus)
library(emmeans)
library(tidyverse)

data <- read_csv("Mariesdata2.csv")

# tidy and wrangel the data to get it in the correct format
tidy_data <- data %>% 
  dplyr::select(Participant, Item, type, 
                Focus, Cleft, YNResp.corr_raw, YNResp.rt_raw, changeResponse.corr_raw) %>%
  filter(type != "filler") %>% 
  dplyr::select(-type)

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

# Results ####
# Viosualisation of all responses
tidy_data %>% 
  group_by(Participant, Focus, Cleft) %>%
  mutate(Response = Response * 100) %>%
  summarise(Response = mean(Response)) %>%
  group_by(Focus, Cleft) %>% 
  summarise(sd = sd(Response), Response = mean(Response)) %>%
  ungroup() %>%
  mutate(Cleft = recode_factor(Cleft, y = "Cleft Present", n = "Cleft Absent")) %>%
  mutate(Focus = recode_factor(Focus, Early = "Early Focus", Late = "Late Focus")) %>%
  ggplot(aes(x = Focus:Cleft, y = Response, fill = Focus:Cleft)) +
  geom_col() +
  geom_errorbar(aes(ymin = Response - sd, ymax = Response + sd), colour="black", width=.1) +
  guides(fill = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Change Detection Rate (%)", title = "All Responses")

# Visualisation of RT data for all responses 
tidy_data %>% 
  mutate(Cleft = recode_factor(Cleft, y = "\nCleft Present", n = "\nCleft Absent")) %>%
  mutate(Focus = recode_factor(Focus, Early = "Early Focus", Late = "Late Focus")) %>%
  ggplot(aes(x = Focus:Cleft, y = RT, colour = Focus:Cleft)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time (ms.)", title = "RT for all Responses")

# Visualisation of RT data when people spot a change (regardless of whether they get it correct)
tidy_data %>% 
  filter(Response == "1") %>%
  mutate(Cleft = recode_factor(Cleft, y = "\nCleft Present", n = "\nCleft Absent")) %>%
  mutate(Focus = recode_factor(Focus, Early = "Early Focus", Late = "Late Focus")) %>%
  ggplot(aes(x = Focus:Cleft, y = RT, colour = Focus:Cleft)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time (ms.)", title = "When a change is spotted (regardless of typed response)")

# Visualisation of RT data when people spot a change AND correctly identify it
tidy_data %>% 
  filter(Correct_word == "1") %>%
  mutate(Cleft = recode_factor(Cleft, y = "\nCleft Present", n = "\nCleft Absent")) %>%
  mutate(Focus = recode_factor(Focus, Early = "Early Focus", Late = "Late Focus")) %>%
  ggplot(aes(x = Focus:Cleft, y = RT, colour = Focus:Cleft)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time (ms.)", title = "Correctly Identified Responses Only")

# Response Rate Analysis ####
# We conducted a generalised mixed model on the change detection rate with Focus 
# and Cleft as fixed effects and participant and item as random intercepts as 
# the fully maximal model would not converge. A summary of the model can be seen
# below. There was a main effect of Focus, but no effect of Cleft, and no 
# interaction between Focus and Cleft.

model <- glmer(Response ~ Focus * Cleft + (1 | Participant) +
                 (1 | Item), data = tidy_data, family = "binomial")
summary(model)
emmeans(model, pairwise ~ Focus, type = "response")

# We conducted a further generalised mixed model on the change detection rate 
# for cases where people spotted a change (regardless of whether they correctly 
# identified it). The fixed and random effects structures are identical to the 
# previous model. A summary of the model can be seen below. There was a main 
# effect of Focus, but no effect of Cleft, and no interaction between Focus and 
# Cleft.

corr_data <- tidy_data %>% filter(Response == "1")
model <- glmer(Correct_word ~ Focus * Cleft + (1 | Participant) +
                 (1 | Item), data = corr_data, family = "binomial")
summary(model)

emmeans(model, pairwise ~ Focus, type = "response")
emmeans(model, pairwise ~ Cleft, type = "response")

# Reaction Time Analysis
# For the reaction time data, three mixed models were built. Due to 
# non-normality of the residuals and convergence issues modeling under the 
# gamma distribution, the mixed models were constructed with the dependent 
# variable being log-transformed reaction times.  The first model included
# reaction times for detected as well as non-detected changes (i.e., all 
# experimental trials regardless of response). This model involved random 
# intercepts for Participants and Items and random slopes for Focus.

# A summary of the model can be seen below. As expected there was an effect of 
# Focus but no effect of Cleft. Unlike the response models, the reaction time 
# model did show an interaction between Focus and Cleft.

all_data <- tidy_data %>% 
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Focus * Cleft + (1 + Focus | Participant) +
                (1 + Focus | Item), data = all_data)
summary(model)

emmeans(model, pairwise ~ Focus * Cleft, adjust = "bonferroni")

# For the second reaction time model, we analysed the data for cases where 
# people spotted a change regardless of wheether they correctly identified it.  
# This time the random effects structure included random intercepts for 
# participants and items, and additive random slopes for our two fixed effects 
# of Focus and Cleft. This again showed an effect of Focus,  no effect of Cleft 
# and an interaction between Focus and Cleft.

# Mixed model on RT data when people spot a change
correct_only_data <- corr_data %>% 
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Focus * Cleft + (1 + Focus + Cleft | Participant) +
                (1 + Focus + Cleft | Item), data = correct_only_data)
summary(model)

emmeans(model, pairwise ~ Focus * Cleft, adjust = "bonferroni")

# Finally, for the third reaction time modelwe analysed the data for cases 
# where people spotted a change AND correctly identified it.  This time the 
# random effects structure included random intercepts for participants and items, 
# and additive random slopes for our two fixed effects of Focus and Cleft for 
# the Participant random effect term but a random slope of Focus for the Item 
# random effect term. This again showed an effect of Focus, no effect of Cleft 
# and an interaction between Focus and Cleft.

# Mixed model on RT data when people spot a change and correctly identify it
correct_only_data <- corr_data %>% 
  filter(Correct_word == "1") %>%
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Focus * Cleft + (1 + Focus + Cleft | Participant) +
                (1 + Focus | Item), data = correct_only_data)
summary(model)

emmeans(model, pairwise ~ Focus * Cleft, adjust = "bonferroni")

# Speed vs. accuracy plot
tidy_data$Response <- as.factor((tidy_data$Response))

tidy_data %>% 
  mutate(Cleft = recode_factor(Cleft, y = "\nCleft Present", n = "\nCleft Absent")) %>%
  mutate(Focus = recode_factor(Focus, Early = "Early Focus", Late = "Late Focus")) %>%
  mutate(Response = recode_factor(Response, "1" = "\nCorrect Response", "0" = "\nIncorrect Response")) %>%
  ggplot(aes(x = Focus:Cleft:Response, y = RT, colour = Response)) +
  geom_jitter(width = .2, alpha = .8, size = .5) +
  guides(colour = FALSE) + 
  scale_colour_hue(h = c(0, 360)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time", title = NULL)
