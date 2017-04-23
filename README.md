# typogenetics

Implementation of Typogenetics from Hofstadter's "Godel, Escher, Bach"

### Simplifying assumptions

- at the start of enzyme application, when there are multiple places that can be bound to, choose the left-most valid binding
- when the ribosomes produce more than one enzyme, apply all of those enzymes to the original strand to produce the set of daughter strands
