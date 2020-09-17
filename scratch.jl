foo=nothing

function x()
  foo=2
  global foo=foo
end
