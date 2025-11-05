import { Suspense } from 'react';
import { Outlet, Route, Routes } from 'react-router-dom';
import DashboardPage from './dashboard/DashboardPage';
import ProjectsPage from './projects/ProjectsPage';
import RunsPage from './runs/RunsPage';
import DataPage from './data/DataPage';
import SettingsPage from './settings/SettingsPage';
import AuthPage from './auth/AuthPage';
import Layout from '../components/Layout';
import LoadingView from '../components/LoadingView';

const App = () => {
  return (
    <Routes>
      <Route element={<AuthPage />} path="/auth" />
      <Route
        element={
          <Layout>
            <Suspense fallback={<LoadingView />}>
              <Outlet />
            </Suspense>
          </Layout>
        }
      >
        <Route element={<DashboardPage />} index />
        <Route element={<ProjectsPage />} path="projects" />
        <Route element={<RunsPage />} path="runs" />
        <Route element={<DataPage />} path="data" />
        <Route element={<SettingsPage />} path="settings" />
      </Route>
    </Routes>
  );
};

export default App;
