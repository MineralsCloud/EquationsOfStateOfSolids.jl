using Preferences: set_preferences!, @load_preference, @set_preferences!

export set_search_interval

get_search_interval() = Tuple(@load_preference("search_interval"))
set_search_interval(search_interval::Tuple{Real,Real}) =
    @set_preferences!("search_interval" => collect(search_interval))

set_preferences!(Inverse, "search_interval" => [eps(), 2]; export_prefs = true)
