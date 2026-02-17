using Documenter, InitialValueProblems

makedocs(
        sitename="InitialValueProblems",
        format=Documenter.HTML(edit_link="main")
)

deploydocs(
      repo = "github.com/vlc1/InitialValueProblems.jl",
      devbranch = "main"
)
