m_ok = function(m)
{
  paste("{crayon::green(clisymbols::symbol$tick)}", m)
}

m_wrn = function(m)
{
  paste("{crayon::yellow('!')}", m)
}

m_err = function(m)
{
  paste("{crayon::red(clisymbols::symbol$cross)}", m)
}

m_inf = function(m)
{
  paste("{crayon::white(clisymbols::symbol$info)}", m)
}

m_txt = function(m, symbol = 'clisymbols::symbol$pointer')
{
  paste('{', symbol, '}', m)
}
