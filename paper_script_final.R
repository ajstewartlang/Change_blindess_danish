library(lme4)
library(lmerTest)
library(fitdistrplus)
library(emmeans)
library(tidyverse)
library(Hmisc)

data <- read_csv("Mariesdata2.csv")

# tidy and wrangle the data to get it in the correct format
tidy_data <- data %>%
  dplyr::select(Participant, Item, type, 
                Focus, Cleft, YNResp.corr_raw, YNResp.rt_raw, changeResponse.corr_raw) %>%
  filter(type != "filler") %>% 
  dplyr::select(-type)

colnames(tidy_data) <- c("Participant", "Item", "Position", "Construction", 
                         "Response", "RT", "Correct_word")

tidy_data[tidy_data$Position == "y",]$Position <- "Early"
tidy_data[tidy_data$Position == "n",]$Position <- "Late"

tidy_data[tidy_data$Construction == "y",]$Construction <- "Cleft"
tidy_data[tidy_data$Construction == "n",]$Construction <- "Non-Cleft"

tidy_data$Position <- as.factor(tidy_data$Position)
tidy_data$Construction <- as.factor(tidy_data$Construction)
tidy_data$Participant <- as.factor(tidy_data$Participant)
tidy_data$Item <- as.factor(tidy_data$Item)

# set up contrasts for the two factors
my_contrasts <- matrix(c(.5, -.5))

contrasts(tidy_data$Position) <- my_contrasts
contrasts(tidy_data$Construction) <- my_contrasts

# convert RT variable to milliseconds
tidy_data$RT <- tidy_data$RT * 1000

# Results ####
# Visualisation of all responses - regardless of correctness
tidy_data %>% 
  group_by(Participant, Position, Construction) %>%
  mutate(Response = Response * 100) %>%
  summarise(Response = mean(Response)) %>%
  group_by(Position, Construction) %>% 
  summarise(sd = sd(Response), Response = mean(Response)) %>%
  ungroup() %>%
  mutate(Position = recode_factor(Position, Early = "Early Position", Late = "Late Position")) %>%
  ggplot(aes(x = Construction:Position, y = Response, fill = Construction:Position)) +
  geom_col() +
  geom_errorbar(aes(ymin = Response - sd, ymax = Response + sd), colour="black", width=.1) +
  guides(fill = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Change Detection Rate (%)", title = "All Responses")

# Visualisation of RT data for all responses - regardless of correctness
tidy_data %>% 
  mutate(Position = recode_factor(Position, Early = "Early Position", Late = "Late Position")) %>%
  ggplot(aes(x = Construction:Position, y = RT, colour = Construction:Position)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time (ms.)", title = "RT for all Responses")

# Visualisation of RT data when people spot a change (regardless of whether they get it correct)
tidy_data %>% 
  filter(Response == "1") %>%
  mutate(Position = recode_factor(Position, Early = "Early Position", Late = "Late Position")) %>%
  ggplot(aes(x = Construction:Position, y = RT, colour = Construction:Position)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time (ms.)", title = "When a change is spotted (regardless of typed response)")

# Visualisation of RT data when people spot a change AND correctly identify it
tidy_data %>% 
  filter(Correct_word == "1") %>%
  mutate(Position = recode_factor(Position, Early = "Early Position", Late = "Late Position")) %>%
  ggplot(aes(x = Construction:Position, y = RT, colour = Construction:Position)) +
  geom_violin() +
  geom_jitter(width = .2, alpha = .2) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  guides(colour = FALSE) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time (ms.)", title = "Correctly Identified Responses Only")

# Response Rate Analysis ####
# We built a generalised mixed model on the change detection rate with Position 
# and Construction as fixed effects and participant and item as random intercepts as 
# the fully maximal model would not converge. A summary of the model can be seen
# below. There was a main effect of Position, but no effect of Construction, and no 
# interaction between Position and Construction.

model <- glmer(Response ~ Position * Construction + (1 | Participant) +
                 (1 | Item), data = tidy_data, family = "binomial")
summary(model)
emmeans(model, pairwise ~ Position, type = "response")

# We conducted a further generalised mixed model on the change detection rate 
# for cases where people spotted a change (regardless of whether they correctly 
# identified it). The fixed and random effects structures are identical to the 
# previous model. There was a main effect of Position, but no effect of Construction, and
# no interaction between Position and Construction.

model <- glmer(Response ~ Position * Construction + (1 | Participant) +
                 (1 | Item), data = tidy_data, family = "binomial")
summary(model)

emmeans(model, pairwise ~ Position, type = "response")

# The next analyisis is when people correctly spot a word has changed, and also
# correctly identify what that change was.
# The fixed and random effects structures are identical to the 
# previous models. A summary of the model can be seen below. There was a main 
# effect of Position, but no effect of Construction, and no interaction between Position and 
# Construction.

corr_data <- tidy_data %>% filter(Response == "1")
model <- glmer(Correct_word ~ Position * Construction + (1 | Participant) +
                 (1 | Item), data = corr_data, family = "binomial")
summary(model)

emmeans(model, pairwise ~ Position, type = "response")
emmeans(model, pairwise ~ Construction, type = "response")

# Reaction Time Analysis
# For the reaction time data, three mixed models were built. Due to 
# non-normality of the residuals and convergence issues modeling under the 
# gamma distribution, the mixed models were constructed with the dependent 
# variable being log-transformed reaction times.  

# The first model included reaction times for detected as well as non-detected
# changes (i.e., all experimental trials regardless of response). This model 
# involved random intercepts and slopes for the factor 'Position' on both Participants 
# and Items.

# A summary of the model can be seen below. As expected there was an effect of 
# Position but no effect of Construction. Unlike the response models, the reaction time 
# model did show an interaction between Position and Construction.

log_data <- tidy_data %>% 
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Position * Construction + (1 + Position | Participant) +
                (1 + Position | Item), data = log_data)
summary(model)

emmeans(model, pairwise ~ Position * Construction, adjust = "bonferroni")

# For the second reaction time model, we analysed the data for cases where 
# people spotted a change regardless of wheether they correctly identified it.  
# This time the random effects structure included random intercepts for 
# participants and items, and additive random slopes for our two fixed effects 
# of Position and Construction. This again showed an effect of Position, no effect
# of Construction and an interaction between Position and Construction.

# Mixed model on RT data when people spot a change
log_spot_data <- tidy_data %>%
  filter(Response == 1) %>%
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Position * Construction + (1 + Position + Construction | Participant) +
                (1 + Position + Construction | Item), data = log_spot_data)
summary(model)

emmeans(model, pairwise ~ Position * Construction, adjust = "bonferroni")

# Finally, for the third reaction time model we analysed the data for cases 
# where people spotted a change AND correctly identified it.  This time the 
# random effects structure included random intercepts for participants and items, 
# and additive random slopes for our two fixed effects of Position and Construction for 
# the Participants and a random slope of Position for Items. 
# This again showed an effect of Position, no effect of Construction 
# and an interaction between Position and Construction.

# Mixed model on RT data when people spot a change and correctly identify it
log_corr_data <- tidy_data %>% 
  filter(Correct_word == "1") %>%
  mutate(logRT = log(RT))

model <- lmer(logRT ~ Position * Construction + (1 + Position + Construction | Participant) +
                (1 + Position | Item), data = log_corr_data)
summary(model)

emmeans(model, pairwise ~ Position * Construction, adjust = "bonferroni")

# Speed vs. accuracy plot
tidy_data$Response <- as.factor((tidy_data$Response))

tidy_data %>% 
  mutate(Construction = recode_factor(Construction, Cleft = "\nCleft", Non-Cleft = "\nNon-Cleft")) %>%
  mutate(Position = recode_factor(Position, Early = "Early Position", Late = "Late Position")) %>%
  mutate(Response = recode_factor(Response, "1" = "\nCorrect Response", "0" = "\nIncorrect Response")) %>%
  ggplot(aes(x = Position:Construction:Response, y = RT, colour = Response)) +
  geom_jitter(width = .2, alpha = .8, size = .5) +
  guides(colour = FALSE) + 
  scale_colour_hue(h = c(0, 360)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "black") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = "Reaction Time", title = NULL)

