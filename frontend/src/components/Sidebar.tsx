const Sidebar = () => {
  return (
    <aside className="sidebar">
      <h1 className="sidebar__logo">Omicsomics</h1>
      <nav className="sidebar__nav">
        <a href="/" className="sidebar__link">
          Dashboard
        </a>
        <a href="/projects" className="sidebar__link">
          Projects
        </a>
        <a href="/runs" className="sidebar__link">
          Runs
        </a>
        <a href="/data" className="sidebar__link">
          Data
        </a>
        <a href="/settings" className="sidebar__link">
          Settings
        </a>
      </nav>
    </aside>
  );
};

export default Sidebar;
