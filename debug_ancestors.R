library(dagitty)
d <- dagitty("dag { Brain -> Migration; Brain -> clutch; Migration -> Lifespan; clutch -> Lifespan }")
print(dagitty::ancestors(d, "Lifespan"))
print(dagitty::descendants(d, "Brain"))
