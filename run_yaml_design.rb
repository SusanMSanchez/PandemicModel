#!/usr/bin/env ruby
require 'yaml'

def is_f?(s)
    !!Float(s) rescue false
end

def is_i?(s)
    /\A[-+]?\d+\z/ === s
end

def convert(s)
  return s.to_i if is_i? s
  return s.to_f if is_f? s
  s
end


design = File.open(ARGV.shift).readlines.map { |line| line.strip.split(/[\s,]+/) }
labels = design.shift.map { |elt| elt.to_sym }

YAML::load_stream( ARGF ) do |model_data|
  design.each do |design_pt|
    labels.each_with_index { |key, i| model_data[key] = convert(design_pt[i]) }
    puts model_data.to_yaml
  end
end
