class TimeVaryingStats

  attr_reader :average

  def initialize
    @last_update_time = 0.0
    @last_update_value = 0.0
    @running_tally = 0.0
  end

  def update(time:, value:)
    elapsed_time = time - @last_update_time
    @running_tally += elapsed_time * @last_update_value
    @last_update_value = value
    @last_update_time = time
    @average = @running_tally / time if time > 0
  end
end

if __FILE__ == $PROGRAM_NAME
  values = [1, 2, 3, 0]
  times = [10, 20, 100, 200]
  stats = TimeVaryingStats.new
  values.each_with_index { |v,i| p stats.update(value: v, time: times[i]) }
  p stats.average
end
