lines = []
len = -1
$stdin.each_line do |l|
    l.strip!
    if len == -1
        len = l.length
    else
        if( len != l.length )
            raise "different line lengths"
        end
    end
    lines << l
end


puts "#{lines.length} #{len}"
i = 0
lines.each do |l|
    puts "#{i} #{l}"
    i+=1
end